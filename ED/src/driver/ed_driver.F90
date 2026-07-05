!==========================================================================================!
!==========================================================================================!
!     Main subroutine that initialises the several structures for the Ecosystem Demography !
! Model 2.  In this version 2.1 all nodes solve some polygons, including the node formerly !
! known as master.                                                                         !
!------------------------------------------------------------------------------------------!
subroutine ed_driver()
   ! ---- Use statements. -------------------------------------------------------------------!
   ! 这些是 use 语句，代表从其他模块（Modules）中导入特定的子程序或变量。
   ! only 关键字是一种安全机制，限制只引入冒号后面的内容，避免命名空间污染。
   ! 例如，从 lsm_hyd（陆面水文模块）中只引入 initHydrology（水文初始化函数）。

   use update_derived_utils , only : update_derived_props          ! ! subroutine
   use lsm_hyd              , only : initHydrology                 ! ! subroutine
   use ed_met_driver        , only : init_met_drivers              & ! subroutine
                                   , read_met_drivers_init         & ! subroutine
                                   , update_met_drivers            ! ! subroutine
   use ed_init_history      , only : resume_from_history           ! ! subroutine
   use ed_init              , only : set_polygon_coordinates       & ! subroutine
                                   , sfcdata_ed                    & ! subroutine
                                   , load_ecosystem_state          & ! subroutine
                                   , read_obstime                  ! ! subroutine
   use grid_coms            , only : ngrids                        & ! intent(in)
                                   , time                          & ! intent(inout)
                                   , timmax                        ! ! intent(inout)
   use ed_state_vars        , only : allocate_edglobals            & ! sub-routine
                                   , filltab_alltypes              & ! sub-routine
                                   , edgrid_g                      ! ! intent(inout)
   use ed_misc_coms         , only : runtype                       & ! intent(in)
                                   , iooutput                      ! ! intent(in)
   use soil_coms            , only : alloc_soilgrid                ! ! sub-routine
   use ed_node_coms         , only : mynum                         & ! intent(in)
                                   , nnodetot                      & ! intent(in)
                                   , sendnum                       ! ! intent(in)
#if defined(RAMS_MPI)
   use ed_node_coms         , only : recvnum                       ! ! intent(in)
#endif
   use detailed_coms        , only : idetailed                     & ! intent(in)
                                   , patch_keep                    ! ! intent(in)
   use phenology_aux        , only : first_phenology               ! ! subroutine
   use hrzshade_utils       , only : init_cci_variables            ! ! subroutine
   use canopy_radiation_coms, only : ihrzrad                       ! ! intent(in)
   use random_utils         , only : init_random_seed              ! ! subroutine
   implicit none ! 强制关闭 Fortran 的隐式类型规则（旧 Fortran 中以 i,j,k,l,m,n 开头的变量默认为整型）
   !----- Included variables. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h' ! MPI commons 当你安装了 MPI 之后，mpif.h 通常会被存放在 MPI 安装目录下的 include 文件夹中
#endif
   !----- Local variables. ----------------------------------------------------------------!
   character(len=12)           :: c0 ! 长度为 12 的字符串，用于格式化输出时间。
   character(len=12)           :: c1 ! 长度为 12 的字符串，用于格式化输出时间。
   integer                     :: ifm ! 循环变量，表示当前处理的网格编号。
   integer                     :: ping ! 用于 MPI 进程间通信的整数变量，通常用于同步或发送信号。
   real                        :: t1 ! 用于记录 CPU 时间。
   real                        :: w1 ! 用于记录开始时间的 wall clock 时间。
   real                        :: w2 ! 用于记录结束时间的 wall clock 时间。
   real                        :: wtime_start ! 用于记录整个初始化过程开始的 wall clock 时间。
   logical                     :: patch_detailed ! 逻辑变量，指示是否只保留一个补丁进行详细分析。
   !----- Local variable (MPI only). ------------------------------------------------------!
#if defined(RAMS_MPI)
   integer                     :: ierr ! 如果是并行编译，声明 ierr 变量，用来接收 MPI 函数执行后的错误状态码（0 代表成功）。
#endif
   !----- External functions. -------------------------------------------------------------!
   real             , external :: walltime    ! wall time 声明 walltime 是一个外部定义的函数（External Function），它返回一个实数，用于精准测量程序运行耗时。
   !---------------------------------------------------------------------------------------!

   ping = 741776 ! 给 ping 赋一个随机的初始整数（用作 MPI 通信的哑变量数据）

   !---------------------------------------------------------------------------------------!
   !      Set the initial time.                                                            !
   !---------------------------------------------------------------------------------------!
   wtime_start = walltime(0.) ! 调用 walltime(0.) 记录当前绝对时间，存入 wtime_start 作为基准。
   w1          = walltime(wtime_start) ! 接着计算当前相对时间 w1。
   !---------------------------------------------------------------------------------------!

   ! 2. 并行随机数种子初始化（Token 环机制）
   !---------------------------------------------------------------------------------------!
   !     Initialise random seed -- the MPI barrier may be unnecessary, added because the   !
   ! jobs may the the system random number generator.                                      !
   !---------------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
   if (mynum /= 1) then ! 并行模式下：如果当前进程号 mynum 不是 1（不是第一个进程），它就会死等
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,79,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   else
      write (unit=*,fmt='(a)') ' [+] Init_random_seed...'
   end if
#else ! 单机模式下：直接打印该行。
      write (unit=*,fmt='(a)') ' [+] Init_random_seed...'
#endif
   call init_random_seed() ! 解析：调用函数初始化随机数发生器。

#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,79,MPI_COMM_WORLD,ierr)
   if (nnodetot /= 1    ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   ! 初始化完随机数后，如果当前进程不是最后一个进程（mynum < nnodetot），就把 ping 信号发送给下一个进程（sendnum），
   ! 解除它的等待状态。最后 MPI_Barrier 让所有进程在这里集合，等所有人都初始化完再一起往下走。
   !---------------------------------------------------------------------------------------!


   ! 3. 生态参数与 XML 配置覆盖
   !---------------------------------------------------------------------------------------!
   !      Set most ED model parameters that do not come from the namelist (ED2IN).         !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ed_Ecosystem_Params...'
   call load_ed_ecosystem_params()
   ! 为了防止多进程同时往屏幕打印造成混乱，这里规定只有最后一个进程（mynum == nnodetot）负责输出提示信息。
   ! 然后所有进程共同执行 load_ed_ecosystem_params() 加载默认生态参数。
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      Overwrite the parameters in case a XML file is provided                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
   if (mynum /= 1) then
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,80,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   else
      write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
   end if
#else
   write (unit=*,fmt='(a)') ' [+] Checking for XML config...'
#endif

   call overwrite_with_xml_config(mynum)
   ! 读取 XML 配置文件并覆盖默认参数，防止多个进程同时读取磁盘引发 I/O 崩溃。

#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,80,MPI_COMM_WORLD,ierr)
   if (nnodetot /= 1 )    call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise any variable that should be initialised after the xml parameters have !
   ! been read.                                                                            !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Init_derived_params_after_xml...'
   call init_derived_params_after_xml()
   ! 计算那些依赖于新参数的派生生态参数
   !---------------------------------------------------------------------------------------!

   !-----Always write out a copy of model parameters in xml--------------------------!
   if (mynum == nnodetot) then 
       write (unit=*,fmt='(a)') ' [+] Write parameters to xml...'      
       call write_ed_xml_config()
       ! 最终生效的参数重新导出一份 XML 备份文件
   endif
   !---------------------------------------------------------------------------------!


   ! 4. 地形、土壤与网格空间初始化
   !---------------------------------------------------------------------------------------!
   !      In case this simulation will use horizontal shading, initialise the landscape    !
   ! arrays.                                                                               !
   !---------------------------------------------------------------------------------------!
   ! 检查 ihrzrad（水平地理辐射标记）。
   !  如果是 0，什么都不做（continue）；
   ! 如果是其他值（说明开启了山地地形遮阴），则初始化地形景观数组（init_cci_variables()）。
   select case (ihrzrad)
   case (0)
      continue
   case default
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Init_cci_variables...'
      call init_cci_variables() ! 初始化地形景观数组（init_cci_variables()
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Allocate soil grid arrays.                                                       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Alloc_Soilgrid...'
   call alloc_soilgrid() ! 分配土壤网格所需的内存空间。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Set some polygon-level basic information, such as lon/lat/soil texture.          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Set_Polygon_Coordinates...'
   call set_polygon_coordinates() ! 设置模拟多边形的基础地理信息（如经纬度）。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialize inherent soil and vegetation properties.                              !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Sfcdata_ED...'
   call sfcdata_ed() ! 加载陆面基础数据（土壤质地、潜在植被覆盖等）。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   if (trim(runtype) == 'HISTORY' ) then
      ! 切掉前后空格后，判断运行类型 runtype 是否为 'HISTORY'（历史断点续跑）。
      !------------------------------------------------------------------------------------!
      !      Initialize the model state as a replicate image of a previous  state.         !
      !------------------------------------------------------------------------------------!
      if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
      if (mynum /= 1) then
         call MPI_RECV(ping,1,MPI_INTEGER,recvnum,81,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      else
         write (unit=*,fmt='(a)') ' [+] Resume_From_History...'
      end if
#else
      write (unit=*,fmt='(a)') ' [+] Resume_From_History...'
#endif
      call resume_from_history()
      ! 如果是续跑，同样使用接力棒机制（Tag=81），确保各节点按顺序读取历史输入文件（resume_from_history()），
      ! 恢复到上一次结束时的状态

#if defined(RAMS_MPI)
      if (mynum < nnodetot ) then
         call MPI_Send(ping,1,MPI_INTEGER,sendnum,81,MPI_COMM_WORLD,ierr)
      end if

      if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      !------------------------------------------------------------------------------------!
   else

      !------------------------------------------------------------------------------------!
      !      Initialize state properties of polygons/sites/patches/cohorts.                !
      !------------------------------------------------------------------------------------!
      if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Load_Ecosystem_State...'
      call load_ecosystem_state() 
      ! 如果是全新运行（else 分支），直接调用 load_ecosystem_state()，根据配置初始化植被、同龄群等状态。
      !------------------------------------------------------------------------------------!
   end if

   !---------------------------------------------------------------------------------------!
   !      In case the runs is going to produce detailed output, we eliminate all patches   !
   ! but the one to be analysed in detail.  Special cases:                                 !
   !  0 -- Keep all patches.                                                               !
   ! -1 -- Keep the one with the highest LAI                                               !
   ! -2 -- Keep the one with the lowest LAI                                                !
   !---------------------------------------------------------------------------------------!
   patch_detailed = ibclr(idetailed,5) > 0
   ! ibclr(idetailed,5) 是 Fortran 的位操作函数，将 idetailed 的第 5 位清零，用其余位判断是否需要输出极度详细的斑块信息。
!   if (patch_detailed) then
      call exterminate_patches_except(patch_keep)
      ! 调用 exterminate_patches_except 杀掉不需要研究的斑块（具体逻辑见第三部分）。
!   end if
   !---------------------------------------------------------------------------------------!

   ! 7. 气象驱动数据初始化
   !---------------------------------------------------------------------------------------!
   !      Initialize meteorological drivers.                                               !
   !---------------------------------------------------------------------------------------!
   ! 解析：又是一套接力棒同步机制（Tag=82），各节点排队初始化并读取气象强迫（强迫场）输入文件
   ! （init_met_drivers, read_met_drivers_init），
   ! 避免多线程抢占磁盘输入流。
#if defined(RAMS_MPI)
   if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   if (mynum == nnodetot-1) sendnum = 0

#if defined(RAMS_MPI)
   if (mynum /= 1) then
      call MPI_RECV(ping,1,MPI_INTEGER,recvnum,82,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
   else
      write (unit=*,fmt='(a)') ' [+] Init_Met_Drivers...'
   end if
#else
   write (unit=*,fmt='(a)') ' [+] Init_Met_Drivers...'
#endif

   call init_met_drivers()
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Read_Met_Drivers_Init...'
   call read_met_drivers_init()


#if defined(RAMS_MPI)
   if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,82,MPI_COMM_WORLD,ierr)
   if (nnodetot /= 1 ) call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise the site-level meteorological forcing.                                !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Update_met_drivers...'
   do ifm=1,ngrids
      call update_met_drivers(edgrid_g(ifm))
   end do
   ! 遍历所有网格（1 到 ngrids），根据当前时间把第一步的气象数据（温、湿、风、辐射等）更新到各个网格结构体 edgrid_g(ifm) 中。
   !---------------------------------------------------------------------------------------!

   ! 8. 水文、物候与最后准备
   !---------------------------------------------------------------------------------------!
   !      Initialize ed fields that depend on the atmosphere.                              !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Ed_Init_Atm...'
   call ed_init_atm()
   ! 初始化受大气状态直接影响的生态系统底层字段。
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !      Initialize hydrology related variables.                                          !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] initHydrology...'
   call initHydrology()
   ! 初始化水文相关变量（如土壤各层初始水分平衡）。
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Initialise some derived variables.  Skip this in case the simulation is resuming !
   ! from HISTORY.                                                                         !
   !---------------------------------------------------------------------------------------!
   if (trim(runtype) /= 'HISTORY' ) then
      do ifm=1,ngrids
         call update_derived_props(edgrid_g(ifm))
      end do
   end if
   ! 如果不是断点续跑，则遍历网格，计算一些由初始状态推导出来的派生属性（续跑时这些属性已经存在于历史文件中，故跳过）。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Initialise drought phenology.  This should be done after the soil moisture has   !
   ! been set up.                                                                          !
   !---------------------------------------------------------------------------------------!
   if (runtype /= 'HISTORY') then
      do ifm=1,ngrids
         call first_phenology(edgrid_g(ifm))
      end do
   end if
   ! 如果不是续跑，结合刚刚算好的土壤水分和气象，计算初期的干旱/季节性物候特征（植被是否展叶或落叶）。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Fill the variable data-tables with all of the state data.  Also calculate the    !
   ! indexing of the vectors to allow for segmented I/O of hyperslabs and referencing of   !
   ! high level hierarchical data types with their parent types.                           !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write (unit=*,fmt='(a)') ' [+] Filltab_Alltypes...'
   call filltab_alltypes
   ! 建立内存结构体和底层数据指针表格的映射，为高维大规模 HDF5 数据的高效分块 I/O 读写做准备。
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Check how the output was configured and determine the averaging frequency.       !
   !---------------------------------------------------------------------------------------!
   if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Find frqsum...'
   call find_frqsum()
   ! 根据用户配置的输出频率和单位，计算出模型内部用于累积平均的时间间隔（frqsum）以及相关的转换因子。
   ! 这个步骤非常重要，因为它决定了模型在时间积分过程中如何处理输出数据的平均和累积。 
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Read obsevation time list if IOOUTPUT is set as non-zero.                        !
   !                                                                                       !
   ! MLO --- Whenever reading ASCII files, it is a good idea to apply MPI barriers, to     !
   !         avoid two nodes accessing the file at the same time (some file systems do not !
   !         like that).                                                                   !
   !---------------------------------------------------------------------------------------!
   if (iooutput /= 0) then
#if defined(RAMS_MPI)
        if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,62,MPI_COMM_WORLD         &
                                     ,MPI_STATUS_IGNORE,ierr)
#endif
        if (mynum == nnodetot) write(unit=*,fmt='(a)') ' [+] Load obstime_list...'
        call read_obstime()
#if defined(RAMS_MPI)
        if (mynum < nnodetot ) call MPI_Send(ping,1,MPI_INTEGER,sendnum,62,MPI_COMM_WORLD  &
                                            ,ierr)
#endif
    end if
   ! 如果用户配置了非零的 iooutput（表示需要按照特定时间点输出数据），则调用 read_obstime() 读取这些时间点列表。
   !---------------------------------------------------------------------------------------!


   ! 10. 打印初始化完成横幅并启动模型
   !---------------------------------------------------------------------------------------!
   !      Get the CPU time and print the banner.                                           !
   !---------------------------------------------------------------------------------------!
   call timing(1,t1)
   w2 = walltime(wtime_start)
   if (mynum == nnodetot) then
      write(c0,'(f12.2)') t1
      write(c1,'(f12.2)') w2-w1
      write(unit=*,fmt='(/,a,/)') ' === Finish initialization; CPU(sec)='//                &
                                  trim(adjustl(c0))//'; Wall(sec)='//trim(adjustl(c1))//   &
                                  '; Time integration starts (ed_master) ==='
   end if
   timing(1,t1) 获取消耗的 CPU 时间。

w2 = walltime(...) 计算流逝的真实世界时间（墙钟时间）。

   ! 将内部数值格式化为字符串，去掉空格后，由最后一个节点向屏幕打印出类似：
   ! === Finish initialization; CPU(sec)=12.50; Wall(sec)=4.20; Time integration starts ... ===
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! STEP 14. Run the model or skip if it is a zero time run.                              !
   !---------------------------------------------------------------------------------------!
   if (time < timmax) then
      call ed_model()
   end if
   ! 如果当前模型时间 time 小于设定的最大模拟截止时间 timmax，正式调用 ed_model() 启动 ED2 的生态时间步长主循环。
   ! 初始化驱动程序结束。
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine finds which frequency the model should use to normalise averaged    !
! variables.  FRQSUM should never exceed one day to avoid build up and overflows.          !
!------------------------------------------------------------------------------------------!
subroutine find_frqsum()
   ! 该程序计算输出变量平均时的归一化分母（frqsum，单位：秒）。
   use ed_misc_coms, only : unitfast        & ! intent(in)
                          , unitstate       & ! intent(in)
                          , isoutput        & ! intent(in)
                          , ifoutput        & ! intent(in)
                          , itoutput        & ! intent(in)
                          , imoutput        & ! intent(in)
                          , iooutput        & ! intent(in)
                          , idoutput        & ! intent(in)
                          , iqoutput        & ! intent(in)
                          , frqstate        & ! intent(in)
                          , frqfast         & ! intent(in)
                          , dtlsm           & ! intent(in)
                          , radfrq          & ! intent(in)
                          , frqsum          & ! intent(out)
                          , frqsumi         & ! intent(out)
                          , dtlsm_o_frqsum  & ! intent(out)
                          , radfrq_o_frqsum ! ! intent(out)
   use consts_coms, only: day_sec

   implicit none 
   !----- Local variables. ----------------------------------------------------------------!
   logical :: fast_output
   logical :: no_fast_output
   ! 声明两个布尔变量，用于标记是否存在快速（高频）数据分析输出。
   !---------------------------------------------------------------------------------------!


   !----- Ancillary logical tests. --------------------------------------------------------!
   fast_output     = ifoutput /= 0 .or. itoutput /= 0 .or. iooutput /= 0
   no_fast_output = .not. fast_output
   ! 如果快速分析输出（ifoutput）、过渡状态输出（itoutput）或特殊点输出（iooutput）中任意一个不为0，
   ! 则 fast_output 为真。no_fast_output 取其反值。
   !---------------------------------------------------------------------------------------!



   if ( no_fast_output .and. isoutput == 0 .and. idoutput == 0 .and. imoutput == 0 .and.   &
        iqoutput == 0  ) then
      write(unit=*,fmt='(a)') '---------------------------------------------------------'
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '  WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! '
      write(unit=*,fmt='(a)') '---------------------------------------------------------'
      write(unit=*,fmt='(a)') ' You are running a simulation that will have no output...'
      frqsum=day_sec ! This avoids the number to get incredibly large.
   ! 如果你把所有的输出开关全部关掉了（没有任何结果输出），程序会疯狂打印一堆 WARNING!，
   ! 并把归一化周期 frqsum 设为一天（day_sec = 86400秒），防止内部变量无限累加导致内存溢出。

   !---------------------------------------------------------------------------------------!
   !    Mean diurnal cycle is on.  Frqfast will be in seconds, so it is likely to be the   !
   ! smallest.  The only exception is if frqstate is more frequent thant frqfast, so we    !
   ! just need to check that too.                                                          !
   !---------------------------------------------------------------------------------------!
   ! 如果开启了日变化均值输出（iqoutput > 0）。若状态输出单位是秒（unitstate == 0），
   ! 则在状态输出频率 frqstate、快速输出频率 frqfast 和一天（86400秒）之间选最小值作为归一化周期；
   ! 否则在 frqfast 和一天之间取最小值。
   elseif (iqoutput > 0) then
      if (unitstate == 0) then
         frqsum = min(min(frqstate,frqfast),day_sec)
      else
         frqsum = min(frqfast,day_sec)
      end if

   !---------------------------------------------------------------------------------------!
   !     Either no instantaneous output was requested, or the user is outputting it at     !
   ! monthly or yearly scale, force it to be one day.                                      !
   !---------------------------------------------------------------------------------------!
   elseif ((isoutput == 0  .and. no_fast_output) .or.                                      &
           (no_fast_output .and. isoutput  > 0 .and. unitstate > 1) .or.                   &
           (isoutput == 0 .and. fast_output .and. unitfast  > 1) .or.                      &
           (isoutput > 0 .and. unitstate > 1 .and. fast_output .and. unitfast > 1)         &
          ) then
      frqsum=day_sec
   ! 如果用户要求的输出是以月、年为单位的长周期输出（单位代码 > 1），或者根本不需要频繁输出，
   ! 强制将变量平均的归一化截止周期锁死为最大值：一天。
   !---------------------------------------------------------------------------------------!
   !    Only restarts, and the unit is in seconds, test which frqsum to use.               !
   !---------------------------------------------------------------------------------------!
   elseif (no_fast_output .and. isoutput > 0) then
      frqsum=min(frqstate,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Only fast analysis, and the unit is in seconds, test which frqsum to use.          !
   !---------------------------------------------------------------------------------------!
   elseif (isoutput == 0 .and. fast_output) then
      frqsum=min(frqfast,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on and both outputs are in seconds or day scales. Choose the minimum      !
   ! between them and one day.                                                             !
   !---------------------------------------------------------------------------------------!
   elseif (unitfast < 2 .and. unitstate < 2) then 
      frqsum=min(min(frqstate,frqfast),day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on but unitstate is in month or years. Choose the minimum between frqfast !
   ! and one day.                                                                          !
   !---------------------------------------------------------------------------------------!
   elseif (unitfast < 2) then 
      frqsum=min(frqfast,day_sec)
   !---------------------------------------------------------------------------------------!
   !    Both are on but unitfast is in month or years. Choose the minimum between frqstate !
   ! and one day.                                                                          !
   !---------------------------------------------------------------------------------------!
   else
      frqsum=min(frqstate,day_sec)
   end if
   ! 其余的分支处理各种秒、天级高频输出共存的情况，核心逻辑就是：永远在用户要求的输出时间频率与“一天”之间找最小值。
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Find some useful conversion factors.                                              !
   ! 1. FRQSUMI         -- inverse of the elapsed time between two analyses (or one day).  !
   !                       This should be used by variables that are fluxes and are solved !
   !                       by RK4, they are holding the integral over the past frqsum      !
   !                       seconds.                                                        !
   ! 2. DTLSM_O_FRQSUM  -- inverse of the number of the main time steps (DTLSM) since      !
   !                       previous analysis.  Only photosynthesis- and decomposition-     !
   !                       related variables, or STATE VARIABLES should use this factor.   !
   !                       Do not use this for energy and water fluxes, CO2 eddy flux, and !
   !                       CO2 storage.                                                    !
   ! 3. RADFRQ_O_FRQSUM -- inverse of the number of radiation time steps since the         !
   !                       previous analysis.  Only radiation-related variables should use !
   !                       this factor.                                                    !
   !---------------------------------------------------------------------------------------!
   frqsumi         = 1.0    / frqsum
   dtlsm_o_frqsum  = dtlsm  * frqsumi
   radfrq_o_frqsum = radfrq * frqsumi
   ! frqsumi：归一化周期的倒数（$1 / frqsum$）。供龙格库塔法（RK4）积分出的物理通量计算时段均值。
   ! dtlsm_o_frqsum：主陆面过程时间步长占归一化周期的比例。供光合、呼吸作用等状态变量平均时使用。
   ! radfrq_o_frqsum：辐射步长占归一化周期的比例。
   !---------------------------------------------------------------------------------------!


   return
end subroutine find_frqsum
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!    This sub-routine eliminates all patches except the one you want to save...  This      !
! shouldn't be used unless you are debugging the code.                                     !
!------------------------------------------------------------------------------------------!
! 第四部分：exterminate_patches_except 逐行解析
! 该程序用于 Debug，在内存中强行抹除除指定斑块以外的所有植物生态斑块。
subroutine exterminate_patches_except(keeppa)
   ! (导入大批代表 ED2 空间多级结构的指针结构体类型)
   use ed_state_vars  , only : edgrid_g           & ! structure
                             , edtype             & ! structure
                             , polygontype        & ! structure
                             , sitetype           & ! structure
                             , patchtype          ! ! structure
   use grid_coms      , only : ngrids             ! ! intent(in)
   use fuse_fiss_utils, only : terminate_patches  ! ! sub-routine
   !----- Arguments -----------------------------------------------------------------------!
   integer                        , intent(in)  :: keeppa
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)                   , pointer     :: cgrid
   type(polygontype)              , pointer     :: cpoly
   type(sitetype)                 , pointer     :: csite
   type(patchtype)                , pointer     :: cpatch
   integer                                      :: ifm
   integer                                      :: ipy
   integer                                      :: isi
   integer                                      :: ipa
   integer                                      :: keepact
   real             , dimension(:), allocatable :: csite_lai
   ! 声明用于遍历 ED2 树状多级空间解算结构的指针
   ! （网格 $\rightarrow$ 多边形 $\rightarrow$ 站点 $\rightarrow$ 斑块）。
   ! csite_lai 是一个一维动态可分配数组，用来临时存放每个斑块的叶面积指数（LAI）。
   !---------------------------------------------------------------------------------------!

   ! 带有标签（如 gridloop）的三层嵌套循环。利用 => 指针赋值符号，逐级将指针指向当前正在处理的网格、多边形和站点。
   gridloop: do ifm=1,ngrids
      cgrid => edgrid_g(ifm)

      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            select case(keeppa)
            case (0)
               ! 检查输入的保留指令：如果是 0，表示保留所有斑块，直接结束子程序退出。
               return
            case (-2)
               !----- Keep the one with the lowest LAI. -----------------------------------!
               allocate(csite_lai(csite%npatches))
               csite_lai(:) = 0.0
               keepm2loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)
                  if (cpatch%ncohorts > 0) csite_lai(ipa) = sum(cpatch%lai)
               end do keepm2loop
               keepact = minloc(csite_lai,dim=1)
               deallocate(csite_lai)
               ! 动态分配一个大小等于当前站点总斑块数（npatches）的数组。
               ! 遍历所有斑块，如果斑块内有植物群落（ncohorts > 0），就把这个斑块所有植物的叶面积指数加总（sum(cpatch%lai)）存进去。
               ! 使用 Fortran 内置函数 minloc 找出整个站点里 LAI 最小的那个斑块的索引，存入 keepact（即决定把它留下来）。随后释放临时数组。
               !---------------------------------------------------------------------------!
            case (-1)
               !----- Keep the one with the lowest LAI. -----------------------------------!
               allocate(csite_lai(csite%npatches))
               csite_lai(:) = 0.0
               keepm1loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)
                  if (cpatch%ncohorts > 0) csite_lai(ipa) = sum(cpatch%lai)
               end do keepm1loop
               keepact = maxloc(csite_lai,dim=1)
               deallocate(csite_lai)
               ! 同上。如果输入是 -1，利用 maxloc 找出整个站点里 LAI 最大的那个斑块留下来。
               !---------------------------------------------------------------------------!
            case default
               !----- Keep a fixed patch number. ------------------------------------------!
               keepact = keeppa
               
               if (keepact > csite%npatches) then
                  write(unit=*,fmt='(a)')       '-----------------------------------------'
                  write(unit=*,fmt='(a,1x,i6)') ' - IPY      = ',ipy
                  write(unit=*,fmt='(a,1x,i6)') ' - ISI      = ',isi
                  write(unit=*,fmt='(a,1x,i6)') ' - NPATCHES = ',csite%npatches
                  write(unit=*,fmt='(a,1x,i6)') ' - KEEPPA   = ',keeppa
                  write(unit=*,fmt='(a)')       '-----------------------------------------'
                  call fail_whale ()
                  call fatal_error('KEEPPA can''t be greater than NPATCHES'                &
                                  ,'exterminate_patches_except','ed_driver.f90')
               end if
            end select

            patchloop: do ipa=1,csite%npatches
               if (ipa == keepact) then
                  csite%area(ipa) = 1.0
               else
                  csite%area(ipa) = 0.0
               end if
            end do patchloop
            ! 遍历该站点内的所有斑块。如果是我们决定留下来的那个斑块（ipa == keepact），
            ! 强行将其占地面积权重设为 1.0（100% 独占）；所有其他斑块的面积直接抹零（0.0）。

            call terminate_patches(csite)
            ! 调用斑块终结函数。该函数在内部会扫描面积为 0 的斑块，释放其内存并从链表中剔除。

         end do siteloop
      end do polyloop
   end do gridloop

   return
end subroutine exterminate_patches_except
!==========================================================================================!
!==========================================================================================!
