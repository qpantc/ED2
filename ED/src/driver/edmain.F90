!==========================================================================================!
!==========================================================================================!
!                    ****** Ecosystem Demography Model -- ED-2.2 ******                    !
!------------------------------------------------------------------------------------------!
!                                                                                          !
! Main program:                                                                            !
!   - Defines execution strategy (sequential or parallel)                                  !
!   - Enrolls processes at defined execution                                               !
!   - Parses command line arguments, looking for namelist file                             !
!   - Dispatches processes according to execution strategy,                                !
!     invoking master/slave processes or full model process                                !
!   - Destroy processes                                                                    !
!                                                                                          !
!------------------------------------------------------------------------------------------!

! main的基本结构
!  ├── 定义变量
!  ├── 设置默认值
!  ├── 读取命令行参数
!  ├── 初始化 MPI
!  ├── 初始化 OpenMP
!  ├── 判断 master/slave
!  ├── 读取 namelist
!  ├── 调用 ed_driver()
!  └── MPI_Finalize()


program main
   !$ use omp_lib
   implicit none ! 禁止隐式变量

   !---------------------------------------------------------------------------------------!
   !      Local constants.                                                                 !
   !---------------------------------------------------------------------------------------!
   !----- Maximum number of input arguments, including MPI own arguments. -----------------!
   integer, parameter                    :: max_input_args       = 63
   !----- Maximum length of each input argument. ------------------------------------------!
   integer, parameter                    :: max_input_arg_length = 256
   !----- Local variables. ----------------------------------------------------------------!
   integer                               :: numarg                  ! actual # input args
   character(len=max_input_arg_length)   :: cargs(0:max_input_args) ! args 
   character(len=2*max_input_arg_length) :: cargx                   ! scratch
   character(len=max_input_arg_length)   :: name_name 
   integer                               :: machsize ! MPI 总进程数：mpirun -np 16 ./ed2 --> machsize = 16
   integer                               :: machnum ! 当前进程编号（rank），从0开始：mpirun -np 16 ./ed2 --> machnum = 0,1,...,15， machnum = 0 的进程是 master 进程
   integer                               :: ipara ! 执行策略：0 - 单进程；1 - master/slave
   integer                               :: icall ! 进程功能：0 - full model；1 - slave on master/slave run
   integer                               :: nslaves ! slave 进程数：如果是 master/slave 模式，nslaves = machsize - 1；如果是单进程模式，nslaves = 0
   integer                               :: isingle ! 是否强制单进程：0 - 正常运行；1 - 不调用 MPI 初始化，强制单进程
   integer                               :: n
   !------ Local variables (MPI only). ----------------------------------------------------!
#if defined(RAMS_MPI)
   integer                               :: ierr ! MPI error code
#endif
   !------ Intrinsic function to return number of arguments (numarg). ---------------------!
   integer                               :: iargc 
   !----- OMP information. ----------------------------------------------------------------!
   integer                               :: max_threads      !<= omp_get_max_threads()
   integer                               :: num_procs        !<= omp_get_num_procs()
   integer                               :: thread
   integer                               :: cpu
   integer, dimension(64)                :: thread_use
   integer, dimension(64)                :: cpu_use
   integer, external                     :: findmycpu
   !------ MPI interface. -----------------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   !      Default settings, which may change depending on the arguments.                   !
   !---------------------------------------------------------------------------------------!
      !----- Namelist. --------------------------------------------------------------------!
      name_name = 'ED2IN' ! 默认的 namelist 文件名，如果用户没有通过命令行参数指定，那么就使用这个默认值
      !----- Process rank and size (default is single process running full model). --------!
      machsize = 0
      machnum  = 0
      !------------------------------------------------------------------------------------!
      ! IPARA: execution strategy; default single process                                  !
      !    0 iff single process (no MPI run or MPI run with a single process)              !
      !    1 iff master-slave processes (MPI run with more than one process)               !
      !------------------------------------------------------------------------------------!
      ipara    = 0
      !------------------------------------------------------------------------------------!
      ! ICALL: this process function on execution strategy; default full model
      !    0 iff full model (no MPI run) or master on MPI run
      !    1 iff slave on MPI run
      !------------------------------------------------------------------------------------!
      icall    = 0

      !------------------------------------------------------------------------------------!
      ! ISINGLE: force non mpi run                                                         !
      !    0 iff normal run with potential MPI                                             !
      !    1 iff no MPI init is called                                                     !
      !                                                                                    !
      ! MLO: Not sure if this will bring conflicts.  I implemented the non-MPI option as   !
      !      a preprocessor feature, and this skips the MPI commands everywhere in the     !
      !      code.  I tried to reconcile this to the best of my knowledge, it would be     !
      !      good if someone who knows MPI better than me could check.                     !
      !------------------------------------------------------------------------------------!
      isingle  = 0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Summary of execution strategy and process function:                                !
      !           ipara=0        ipara=1                                                   !
      ! icall=0   full model     master on master/slave run                                !
      ! icall=1   impossible     slave  on master/slave run                                !
      !------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Get input arguments.                                                              !
   !---------------------------------------------------------------------------------------!
   numarg=iargc() ! 获取命令行参数个数，不包括程序名称本身：mpirun -np 16 ./ed2 -f ED2IN --> numarg = 3，分别是 "-f", "ED2IN" 和 "-np 16"
   if (numarg > max_input_args) then
      write(unit=*,fmt='(a)')       '-----------------------------------------------------'
      write(unit=*,fmt='(a,1x,i6)') ' NUMARG         = ',numarg
      write(unit=*,fmt='(a,1x,i6)') ' MAX_INPUT_ARGS = ',max_input_args
      write(unit=*,fmt='(a)')       '-----------------------------------------------------'
      call fatal_error('Input argument list length (NUMARG) exceeds MAX_INPUT_ARGS!'       &
                      ,'main','edmain.F90')
   end if
   do n=0,numarg
      call ugetarg(n,cargx) ! 获取第 n 个命令行参数，存储在 cargx 中：mpirun -np 16 ./ed2 -f ED2IN --> n=0 时 cargx = "-f"，n=1 时 cargx = "ED2IN"，n=2 时 cargx = "-np 16"
      if (len_trim(cargx) > max_input_arg_length) then
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         write(unit=*,fmt='(a,1x,i6)') ' ARGUMENT NUMBER      = ',n
         write(unit=*,fmt='(a,1x,a)' ) ' ARGUMENT             = ',cargx
         write(unit=*,fmt='(a,1x,i6)') ' LEN_TRIM(ARGUMENT)   = ',len_trim(cargx)
         write(unit=*,fmt='(a,1x,i6)') ' MAX_INPUT_ARG_LENGTH = ',max_input_arg_length
         write(unit=*,fmt='(a)')       '--------------------------------------------------'
         call fatal_error('Input argument length exceeds MAX_INPUT_ARG_LENGTH!'            &
                         ,'main','edmain.F90')
      end if
      cargs(n)=trim(cargx)//char(0) ! 存储参数到 cargs 数组中，去掉前后空格，并在末尾添加一个 null 字符，
      !！ //char(0) 是 Fortran 中字符串连接的语法，char(0) 表示 ASCII 码为 0 的字符，即 null 字符
      !！ 这样做的目的是为了方便后续使用 C 风格的字符串处理函数，因为 C 语言中的字符串是以 null 字符结尾的，而 Fortran 中的字符串则不一定是 null 结尾的。
   end do
   numarg=numarg+1
   !---------------------------------------------------------------------------------------!

! 这是“编译前预处理指令”，不是 Fortran 运行语句。
! 加#是为编译准备的
#if defined(RAMS_MPI) ! 如果定义了 RAMS_MPI 预处理器宏，说明这是一个 MPI 版本的编译，那么就执行以下代码块
   ! 即当 使用了 -D RAMS_MPI 编译选项时，才会执行以下代码块，否则整个 MPI 代码不会编译。
   !---------------------------------------------------------------------------------------!
   !      Find out if sequential or MPI run; if MPI run, enroll this process.  If          !
   ! sequential execution, machnum and machsize return untouched (both zero); if MPI       !
   ! execution, machnum returns process rank and machsize process size.                    !
   !---------------------------------------------------------------------------------------!
   do n = 1, numarg
       select case (cargs(n)(1:2))
       case ('-s') ! 如果命令行参数中包含 "-s"，说明用户想要强制单进程运行，那么就即使即使编译了 MPI，也不调用 MPI
          isingle = 1
       case default
          continue
       end select
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Decide whether to call MPI.                                                       !
   !---------------------------------------------------------------------------------------!
   select case (isingle) ! 多分支条件判断结构，根据 isingle 的值来决定是否调用 MPI 初始化
   case (0)
      call MPI_Init(ierr) ! 启动 MPI
      call MPI_Comm_rank(MPI_COMM_WORLD,machnum,ierr) ! 获取当前进程的 rank（编号），存储在 machnum 中
      call MPI_Comm_size(MPI_COMM_WORLD,machsize,ierr) ! 获取总的进程数，存储在 machsize 中
   case default
      machnum  = 0
      machsize = 1
   end select
   !---------------------------------------------------------------------------------------!
#else
   !---------------------------------------------------------------------------------------!
   !   Set dummy values for all MPI variables.                                             !
   !---------------------------------------------------------------------------------------!
   isingle       = 1
   machnum       = 0
   machsize      = 1
   num_procs     = 1
   !---------------------------------------------------------------------------------------!
#endif
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Check OMP thread and processor use and availability.  This is done outside the    !
   !  preprocessor "if block" because it is possible to activate multithread without       !
   ! compiling with mpif90.                                                                       !
   !                                                                                       !
   ! Note: One could use omp_get_num_threads() in loop, but that would depend on how many  !
   ! threads were open at the time of its call.                                            !
   !---------------------------------------------------------------------------------------!
   max_threads   = 1
   num_procs     = 1
   thread        = 1
   cpu           = 1
   thread_use(:) = 0
   cpu_use(:)    = 0
   !$ max_threads = omp_get_max_threads() ！ 获取 OpenMP 线程的最大数量，存储在 max_threads 中
   !$ num_procs   = omp_get_num_procs() ！ 获取可用的处理器数量，存储在 num_procs 中

   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(thread,cpu)
   ! OpenMP 并行区域，使用 PARALLEL DO 指令来并行执行循环，
   ! DEFAULT(SHARED) 表示默认情况下变量是共享的，
   ! PRIVATE(thread,cpu) 表示 thread 和 cpu 变量在每个线程中都是私有的
   do n = 1,max_threads
     !$ thread = omp_get_thread_num() + 1
     !$ cpu    = findmycpu() + 1
     thread_use(thread) = 1 ! 标记当前线程被使用了
     cpu_use(cpu)       = 1 ! 标记当前 CPU 被使用了
   end do
   !$OMP END PARALLEL DO
   !---------------------------------------------------------------------------------------!



   !------ Always print the banner. -------------------------------------------------------!
   ! 输出并行信息
   write (*,'(a)')       '+---------------- MPI parallel info: --------------------+'
   write (*,'(a,1x,i6)') '+  - Machnum  =',machnum
   write (*,'(a,1x,i6)') '+  - Machsize =',machsize
   write (*,'(a)')       '+---------------- OMP parallel info: --------------------+'
   write (*,'(a,1x,i6)') '+  - thread  use: ', sum(thread_use)
   write (*,'(a,1x,i6)') '+  - threads max: ', max_threads
   write (*,'(a,1x,i6)') '+  - cpu     use: ', sum(cpu_use)
   write (*,'(a,1x,i6)') '+  - cpus    max: ', num_procs
   write (*,'(a)')       '+  Note: Max vals are for node, not sockets.'
   write (*,'(a)')       '+--------------------------------------------------------+'
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    If this is MPI run, define master or slave process, otherwise, keep default        !
   ! (single process does full model).                                                     !
   !---------------------------------------------------------------------------------------!
   select case (machsize) 
   ! 根据 machsize 的值来判断是单进程运行还是 MPI 运行，
   ! 如果 machsize > 1，说明是 MPI 运行，那么就根据 machnum 的值来判断是 master 进程还是 slave 进程
   case (2:)
      ipara = 1 ! 1 表示 master/slave 模式
      select case (machnum)
      case (0) ! master 进程
         icall = 0
      case default ! slave 进程
         icall = 1
      end select
   case default
      ipara = 0 ! 0 表示单进程模式
      icall = 0
   end select
   !---------------------------------------------------------------------------------------!


   !----- Master process gets number of slaves and sets process ID. -----------------------!
   select case (icall)
   ! 如果 icall = 0，说明这是 master 进程，那么就计算 slave 进程的数量 
   ! nslaves = machsize - 1，因为 master 进程不算在内
   case (0)
      nslaves=machsize-1
   case default
      continue
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Master process parse command line arguments looking for "-f <namelist filename>". !
   !---------------------------------------------------------------------------------------!
   select case (icall)
   ! 主进程解析命令行参数，寻找 "-f <namelist filename>" 这样的参数对，
   ! 如果 icall = 0，说明这是 master 进程，那么就解析命令行参数，
   ! 寻找 "-f <namelist filename>" 这样的参数对，
   case (0)
      do n = 1, numarg
         select case (cargs(n)(1:2))
         case ('-f')
            name_name = cargs(n+1)(1:len_trim(cargs(n+1))-1)
         end select
      end do
   case default
      continue
   end select
   !---------------------------------------------------------------------------------------!



   !----- Read the namelist and initialize the variables in the nodes if needed. ----------!
   select case (icall)
   ! 如果 icall = 0，说明这是 master 进程，那么就调用 ed_1st_master() 函数来读取 namelist 文件并初始化变量，
   ! 这个函数会在 master 进程上执行，并且会把必要的信息传递给 slave 进程，让它们也能正确地执行后续的计算
   case (0)
      call ed_1st_master(ipara,machsize,nslaves,machnum,max_threads,name_name)
   case default ! 如果 icall = 1，说明这是 slave 进程，那么就调用 ed_1st_node() 函数来等待 master 分配任务，
      call ed_1st_node() ! 等待 master 分配任务。
   end select
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Call the main driver: it allocates the structures, initializes the variables,     !
   ! calls the timestep driver, deals with I/O, cooks, does the laundry and irons.  The    !
   ! stand-alone driver tells the master node that it actually has to get a job and do     !
   ! something with its life.  So the driver is passed a zero here, which tells the node   !
   ! with mynum = nnodetot-1 that the master is next in line for sequencing.               !
   !---------------------------------------------------------------------------------------!
   ! 真正计算生态过程的主函数，这里才开始计算：
   ! timestep loop
   ! photosynthesis
   ! respiration
   ! hydrology
   ! demography
   ! disturbance
   ! carbon cycle


   call ed_driver() 


   
   ! 前面没有use ed_driver模块，但是这里能直接调用 ed_driver() 函数，是因为 ed_driver() 函数在 ed_driver.F90 文件中定义了，并且 ed_driver.F90 文件被包含在了编译过程中，所以编译器能够找到这个函数的定义。
   ! 默认能找到 ed_driver() 函数，因为它在 ed_driver.F90 文件中定义了，这个函数会根据 master/slave 模式来执行不同的计算任务，master 进程负责协调和分配任务，而 slave 进程则负责执行具体的计算
   !---------------------------------------------------------------------------------------!



   !----- Finishes execution. -------------------------------------------------------------!
#if defined(RAMS_MPI)
   select case (isingle)
   case (0)
      call MPI_Finalize(ierr) ! 结束 MPI
   case default !都不匹配时
      continue
   end select
#endif
   select case (icall)
   case (0)
      write(unit=*,fmt='(a)') ' ------ ED-2.2 execution ends ------'
   case default
      continue
   end select
   !---------------------------------------------------------------------------------------!

   stop
end program main
!==========================================================================================!
!==========================================================================================!
