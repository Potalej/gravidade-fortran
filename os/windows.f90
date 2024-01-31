module os
  IMPLICIT NONE
  PRIVATE
  PUBLIC listdir
contains

  SUBROUTINE listdir (cmd_ls)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: cmd_ls
    cmd_ls = "dir"    
  END SUBROUTINE listdir

end module os