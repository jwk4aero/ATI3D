#ifdef SS_DOUBLE_PREC
#define SS_PREC dp
#else ! not SS_DOUBLE_PREC
#define SS_PREC sp
#endif ! not SS_DOUBLE_PREC


 module SSVAR
   use stdtypes
!cpb   -------------- SUBSPACES -----------------------------------------
   logical ss_flag
   integer, allocatable,save :: ss_opt(:)
   integer, allocatable,save :: ss_capt(:)
   character*100, allocatable,save :: ss_tag(:)
   character*100, allocatable,save :: ss_tag_grid(:)
   logical, allocatable,save :: ss_blocks(:)
   integer, allocatable,save :: ss_num_var(:)
   character*12,parameter:: ss_path='./Subspaces/'
   integer, allocatable,save :: ss_dim_block(:,:)
   integer, allocatable,save :: ss_dim_proc(:,:)
   integer, allocatable,save :: ss_dim_spec(:,:)
   integer, allocatable,save :: ss_spec_block(:,:,:)
   integer ss_num,ss_num_var_max
   logical, allocatable,save :: ss_in_subspace(:)
   integer, allocatable,save :: ss_communicator(:)
   integer, allocatable,save :: ss_rank(:)
   integer, allocatable :: ss_status(:,:),ss_fh(:)
   integer, allocatable :: ss_mb_specs(:,:)
   integer, allocatable :: ss_blockid(:)
   integer, allocatable :: ss_wall_opt(:)
#ifdef SS_MULTIBLOCK
   integer, allocatable,save :: ss_mb_offset(:)
#endif ! SS_MULTIBLOCK
#ifndef SS_BLOCKING_WRITE
   real(kind=SS_PREC),allocatable :: ss_qio(:)
   logical,allocatable :: ss_first_write(:)
#endif
!cpb   -------------- END   SUBSPACES ----------------------------------
 end module SSVAR


