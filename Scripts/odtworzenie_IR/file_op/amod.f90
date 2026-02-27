      
      module cv
      integer :: rec_lenght
      integer :: position,pos1,pos2,pos3,pos4,num_tokens,ntemp,row_min=0
      logical(1), allocatable :: rows(:),columns(:),mask_1d(:),mask_2d(:,:)
      logical :: row_columns_corrections,l1,l2,center_data,diagonalize_correlation_matrix
      character(20000), pointer :: inp_file(:)
      character(200), pointer :: conf_file(:)
      character(150), pointer :: rows_columns(:)
      character(150), allocatable :: string_tokenized(:)
      character(150) :: atmp,btmp,ctmp,form,inp_file_name
      character(20), allocatable :: sample_table(:),sample_table_red(:)
      integer :: react_coord(2), prod_coord(2), ts_coord,nlines,n_series_red,n_measurements_red,ier
      real(8), allocatable :: r8tmp(:,:),all_spectra_value(:,:),aves_column(:),aves_row(:), &
      spectra_value_ck(:,:),spectra_value_cw(:,:),selected_spectra_value(:,:),sdt_dev_col(:),sdt_dev_row(:),cov_matrix(:,:), &
      correlation_matrix(:,:),x(:),y(:),eigvalues_cov(:),eigen_vec_cov(:,:), &
      xtmp_2d(:,:),s(:,:),b(:,:),b_matrix(:,:),factors_nc(:,:),factors_cc(:,:),factors_rc(:,:),scores(:,:),loadings(:,:), &
      factor_loadings(:,:), factor_scores(:,:),residuals(:,:),singul(:),vi(:),trv(:),erv(:),spectra_odtworzone(:,:)
      real(8) :: skw_xauto,singul0,dfij
!
!  spectra_value_ck - centrowane kolumnowo
!               _cw - centrowane wierszowo
!
      character(20), allocatable :: series_names(:),series_names_red(:)
!
!  NIPALS
!
      real(8), parameter :: eps=1.d-13
      integer, parameter :: max_iter=10000
!
      interface get_file
        subroutine closed_file2table(filename,ndim,table,rec_lenght)
        integer, intent(in) :: ndim,rec_lenght
        character(rec_lenght), pointer :: table(:)
        character(*), intent(in) :: filename
        end subroutine closed_file2table
!
        subroutine opened_file2table(num_unit,ndim,table,rec_lenght)
        integer, intent(in) :: ndim,rec_lenght
        character(rec_lenght), pointer :: table(:)
        end subroutine opened_file2table
      end interface get_file
!
      interface read_config_line
        subroutine read_rc(string,real_value)
        character(*), intent(in) :: string
        real, intent(out) :: real_value
        end subroutine read_rc

        subroutine read_ic(string,int_value)
        character(*), intent(in) :: string
        integer, intent(out) :: int_value
        end subroutine read_ic

        subroutine read_lc(string,logical_value)
        character(*), intent(in) :: string
        logical, intent(out) :: logical_value
        end subroutine read_lc

        subroutine read_charc(string,char_value)
        character(*), intent(in) :: string,char_value
        end subroutine read_charc
      end interface read_config_line
      interface ave_table
        function average_table(table)
        real(8), intent(in) :: table(:)
        end function average_table
      end interface ave_table
      end module cv
!
!*****************************************************************************
!
      module upplow
      public; save
      character(1),dimension(26) :: upper,lower
      character(1),dimension(22) :: no_num
      data upper/'A','B','C','D','E','F','G','H','I','J','K','L','M', &
      'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
      data lower/'a','b','c','d','e','f','g','h','i','j','k','l','m',  &
      'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data no_num /',','/','\','[',']','|','?','>','!','@','#','$','%','^','&','*', &
      '-','_','=','+','{','}'/
      character(1),dimension(17) :: in_numbers
      data in_numbers/'1','2','3','4','5','6','7','8','9','0','e','E','d','D','.',  &
       '-','+'/
      end module upplow
!
!*****************************************************************************
!
      module parse
      integer, allocatable :: sub_sections(:,:)
    end module parse

module const
  ! SP: (4), DP: (8)
  integer,     parameter :: SP = kind(1.0)
  integer(SP), parameter :: DP = selected_real_kind(2 * precision(1.0_SP))
end module const

module comp
  use const
  implicit none
  private
  public :: calc_r,calc_cov

contains
  ! 
  !
  ! :param(in) real(8) x(:): X 
  ! :param(in) real(8) y(:): Y 
  ! :return    real(8)    r: 
  function calc_r(x, y) result(r)
    implicit none
    real(DP), intent(in) :: x(:), y(:)
    real(DP)    :: r
    integer(SP) :: size_x, size_y, i
    real(DP)    :: mean_x, mean_y, cov, var_x, var_y

    size_x = size(x)
    size_y = size(y)
    if (size_x == 0 .or. size_y == 0) then
      print *, "[ERROR] array size == 0"
      stop
    end if
    if (size_x /= size_y) then
      print *, "[ERROR] size(X) != size(Y)"
      stop
    end if

    r = 0.0_DP
    mean_x = sum(x) / size_x
    mean_y = sum(y) / size_y
    cov = sum((x(1:size_x) - mean_x) * (y(1:size_x) - mean_y))
    var_x = sum((x(1:size_x) - mean_x) * (x(1:size_x) - mean_x))
    var_y = sum((y(1:size_x) - mean_y) * (y(1:size_x) - mean_y))
    r = (cov / sqrt(var_x)) / sqrt(var_y)
  end function calc_r
!
  function calc_cov(x, y) result(r)
    implicit none
    real(DP), intent(in) :: x(:), y(:)
    real(DP)    :: r
    integer(SP) :: size_x, size_y, i
    real(DP)    :: mean_x, mean_y, cov, var_x, var_y

    size_x = size(x)
    size_y = size(y)
    if (size_x == 0 .or. size_y == 0) then
      print *, "[ERROR] array size == 0"
      stop
    end if
    if (size_x /= size_y) then
      print *, "[ERROR] size(X) != size(Y)"
      stop
    end if

    r = 0.0_DP
    mean_x = sum(x) / size_x
    mean_y = sum(y) / size_y
    cov = sum((x(1:size_x) - mean_x) * (y(1:size_x) - mean_y))/(size_x-1)
    var_x = sum((x(1:size_x) - mean_x) * (x(1:size_x) - mean_x))
    var_y = sum((y(1:size_x) - mean_y) * (y(1:size_x) - mean_y))
    !r = (cov / sqrt(var_x)) / sqrt(var_y)
    r=cov
  end function calc_cov
end module comp