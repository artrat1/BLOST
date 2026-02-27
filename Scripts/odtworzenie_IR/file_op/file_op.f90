
!  file_op.f90 
!
!  ALL_SPECTRA_VALUE - wszystkie dane z pliku Zeszyt1.csv, procz nazw serii (pierwszy wiersz) 
!  i czestotliwosci (kolumna A)
!
!  SELECTED_SPECTRA_VALUE - czêœæ lub ca³osc danych z pliku ALL_SPECTRA_VALUE. To, ktore wiersze i kolumny
!  sa uwzgledniane opisane jest w pliku INP.DAT. Je¿eli brak pliku INP.DAT lub w odpowiednich polach tego
!  pliku s¹ zera, wtedy SELECTED_SPECTRA_VALUE jest iudentyczny z ALL_SPECTRA_VALUE
!
!  n_measurements (n_measurements_red) (liczby) - iloœæ czêstotliwoœci (wiersze!) pe³na (zredukowana)
!  n_series (n_series_red) (liczby) - iloœæ serii pomiarowych (kolumny!!) pe³na (zredukowana)
!  SAMPLE_TABLE, SAMPLE_TABLE_RED - pe³ny (zredukowany) zbiór próbek (tablice, etykiety wierszy)
!  spectra_value_cw - wek  
!

    program file_op
    use cv; use comp
    implicit real(8) (a-h,o-z)
    logical :: exist
!
!  Input file name and other data from CONFIG file
!
    inp_file_name='Zeszyt1.csv'
    inquire(file='config', exist=exist)
    if(exist)then
        rec_lenght=len(conf_file(1))
        call get_file('config',nlines,conf_file,rec_lenght)
        call process_config_file
      else
        center_data=.true.
    end if
!
100 nlines=0
    rec_lenght=len(inp_file(1))
    call get_file(trim(inp_file_name),nlines,inp_file,rec_lenght)
!
! In the .csv files numbers are separated with semicolons
! Remove these semicolons first
!
    do i=1,nlines
      ll=len_trim(inp_file(i))
      do j=1,ll
         if(inp_file(i)(j:j)==';')inp_file(i)(j:j)=' '
      end do
    end do
!
!  Then allocate working arrays
!
    if(allocated(sample_table))deallocate(sample_table)
    allocate(sample_table(nlines-1)); sample_table(:)(:)=' '
!
!  Find a number of series (tokenize 1-st row)
!    
    call tokenize_string(inp_file(1))
    n_series=num_tokens-1
    n_measurements=nlines-1
    if(allocated(series_names))deallocate(series_names)
    allocate(series_names(n_series))
    series_names=string_tokenized(2:)
    if(allocated(all_spectra_value))deallocate(all_spectra_value)
    allocate(all_spectra_value(n_measurements,n_series)); all_spectra_value=0.d00
!
!   Now read intensities, frequency first (1-st column)
!
    do i=2,nlines
        !print*, i
        read(inp_file(i),*,end=2,err=2)sample_table(i-1),(all_spectra_value(i-1,j),j=1,n_series)
    end do
2    continue
!
!  Which rows and colums are considered
!
    if(allocated(rows))deallocate(rows); allocate(rows(n_measurements)); rows=.true.
    if(allocated(columns))deallocate(columns); allocate(columns(n_series)); columns=.true.
    if(allocated(selected_spectra_value))deallocate(selected_spectra_value)
!
!  Create SELECTED_SPECTRA_VALUE table - only data specified in inp.dat file
!
    call process_inp
    SELECTED_DATA_ONLY:if(row_columns_corrections)then
        n_measurements_red=count(rows)
        n_series_red=count(columns)
        allocate(series_names_red(n_series_red));series_names_red(:)=' '
        allocate(selected_spectra_value(n_measurements_red,n_series_red));selected_spectra_value=0.d00
!
! Columns (series names) first
!
        if(allocated(series_names_red))deallocate(series_names_red);allocate(series_names_red(n_series_red))
        series_names_red(:)(:)=' '
        series_names_red=pack(series_names,mask=columns)
!
! Now rows (measurements) 
!
        if(allocated(sample_table_red))deallocate(sample_table_red);allocate(sample_table_red(n_measurements_red))
        sample_table_red(:)(:)=' '
        sample_table_red=pack(sample_table,mask=rows)
!
!  Now spectra table
!
       if(allocated(mask_2d))deallocate(mask_2d);allocate(mask_2d(n_measurements,n_series))  
       mask_2d=.false.
       FORALL (i=1:n_measurements,j=1:n_series) mask_2d(i,j) = rows(i) .and. columns(j)
       selected_spectra_value=reshape(pack(all_spectra_value,mask=mask_2d),(/n_measurements_red,n_series_red/))
    else
        n_measurements_red=n_measurements
        n_series_red=n_series
        allocate(selected_spectra_value(n_measurements_red,n_series_red))
        selected_spectra_value=all_spectra_value
    end if SELECTED_DATA_ONLY
!
!  Write reduced spectra data to disk file SELECTED_SPECTRA.DAT
!
    open(1,file='selected_spectra.dat')
    form(:)=' '
    form(1:1)='('
    write(form(2:),'(i3)')n_series_red
    form = trim(form) // '(1x,a9))'
    write(1,trim(form))series_names_red
!
    form(:)=' '
    form(1:5)='(a,'
    write(form(7:),'(i3)')n_series_red
    form = trim(form) // '(1x,f9.4))'
    do i=1,n_measurements_red
        write(1,form)trim(sample_table_red(i)),selected_spectra_value(i,1:n_series_red)
    end do 
    close(1,status='keep')
!
!  Center data
!
    call center
!
!  Correlation matrix; X,Y - work arrays
!
    if(allocated(cov_matrix))deallocate(cov_matrix); allocate(cov_matrix(n_series_red,n_series_red))
    if(allocated(correlation_matrix))deallocate(correlation_matrix); allocate(correlation_matrix(n_series_red,n_series_red))
    cov_matrix=1.d00; correlation_matrix=1.d00
    if(allocated(x))deallocate(x,y)
    allocate(x(n_measurements_red),y(n_measurements_red)); x=0.d00; y=0.d00
    do i=1,n_series_red
        do j=1,n_series_red
           !if(i==j)cycle
           !x=selected_spectra_value(:,i)
           !y=selected_spectra_value(:,j)
           x=spectra_value_ck(:,i)
           y=spectra_value_ck(:,j)
          ! cov_matrix(i,j)=calc_r(x, y)
           cov_matrix(i,j)=calc_cov(x, y)
           correlation_matrix(i,j)=calc_r(x, y)
        end do
    end do
   ! cov_matrix=matmul(
    open(1,file='pca_data.dat')
    write(1,*)
    write(1,*)'Data read from file: ',trim(inp_file_name)
    write(1,*)
    if(center_data)then
        write(1,*)'Data are centered'
     else
        write(1,*)'Data are NOT centered'
     end if  
    write(1,*)
    if(diagonalize_correlation_matrix)then
        write(1,*)'CORRELATION matrix is diagonalized'
     else
        write(1,*)'COVARIANCE matrix is diagonalized'
    end if 
    write(1,33)
33  format(/80('-')/)
    write(1,*)'PCA parameters'
    write(1,*)
    write(1,*)'Covariance matrix'
    do i=1,n_series_red
       write(1,'(10000f9.3)')(cov_matrix(i,j),j=1,n_series_red)
    end do
    write(1,*)
    write(1,*)'Correlation matrix'
    do i=1,n_series_red
       write(1,'(10000f9.3)')(correlation_matrix(i,j),j=1,n_series_red)
    end do
    
!
!  Eigenvalues and eigenvectors
!
    call alloc_eigrs
! (A,N,JOBN,D,Z,IZ,WK,IER)
    jobn=12;n=n_series_red
    if(diagonalize_correlation_matrix)then
       cov_matrix=correlation_matrix
    end if
    call EIGRS(cov_matrix,n,jobn,eigvalues_cov,eigen_vec_cov,n_series_red,x,ier)
   ! eigen_vec_cov=cov_matrix
    abserr=1.d-15
   ! call Jacobi(cov_matrix,eigen_vec_cov,abserr,n)
   ! do i=1,n
     ! eigvalues_cov(i)=cov_matrix(i,i)
    !end do
    write(1,*)
    write(1,3)n_series_red,(eigvalues_cov(i),i=1,n_series_red)
3   format('N= ',i5,' Eigenvalues bef. sorting: ',10000f9.3)
    write(1,*)
    write(1,*)'Eigenvectors bef. sorting and resigning: '
    do i=1,n_series_red
       write(1,'(10000f9.3)')(eigen_vec_cov(i,j),j=1,n_series_red)
    end do
!
!  sort eigenvalues and eigenvectors in descending order, EIGRS returns them in ascending order
!
    eigvalues_cov(n_series_red:1:-1)=eigvalues_cov(1:n_series_red)
    eigen_vec_cov(:,n_series_red:1:-1)=eigen_vec_cov(:,1:n_series_red)
!
!  Change sign of an eigenvectors for which more than half components ane  <0
!
    do i=1,n_series_red
       n_positive=count(eigen_vec_cov(:,i)>0.d00)
       n_negative=count(eigen_vec_cov(:,i)<=0.d00)
       if(n_positive>=n_negative)cycle
       eigen_vec_cov(:,i)=eigen_vec_cov(:,i)*(-1.d00)
    end do
    write(1,*)
    write(1,4)n_series_red,(eigvalues_cov(i),i=1,n_series_red)
4   format('N= ',i5,' Eigenvalues after sorting: ',10000f9.3)
    write(1,*)
    write(1,*)'Eigenvectors after sorting and resigning: '
    do i=1,n_series_red
       write(1,'(10000f9.3)')(eigen_vec_cov(i,j),j=1,n_series_red)
    end do
!
!  Matrix L^(1/2)
!
    do i=1,n_series_red
       xtmp_2d(i,i)=dsqrt(eigvalues_cov(i))
    end do
    s=matmul(eigen_vec_cov,xtmp_2d)
    write(1,*)
    write(1,*)'Matrix S^(-1/2): '
    do i=1,n_series_red
       write(1,'(10000f9.3)')(s(i,j),j=1,n_series_red)
    end do
!
!  Matrix B=EV*L(-1/2)
!
    xtmp_2d=0.d00
    do i=1,n_series_red
       xtmp_2d(i,i)=1.d00/dsqrt(eigvalues_cov(i))
    end do
    b_matrix=matmul(eigen_vec_cov,xtmp_2d)
    write(1,*)
    write(1,*)'Matrix B=EV*L(-1/2): '
    do i=1,n_series_red
       write(1,'(10000f9.3)')(b_matrix(i,j),j=1,n_series_red)
    end do
!
! The factor scores from Z*B, where Z is X converted to standard score form.
!
    factors_nc=matmul(selected_spectra_value,b_matrix)
   ! write(1,*)
   ! write(1,*)'Factor scores (NON centered data): '
   ! do i=1,n_measurements_red
   !    write(1,'(10000f9.3)')(factors_nc(i,j),j=1,n_series_red)
   ! end do
!
! Z Matrix - normalized data
!
    write(1,*)
    write(1,*)'Z matrix (COLUMN centered data): '
    do i=1,n_measurements_red
       write(1,'(10000f9.3)')(spectra_value_ck(i,j),j=1,n_series_red)
    end do   
!
    factors_cc=matmul(spectra_value_ck,b_matrix)
    write(1,*)
    write(1,*)'Factor scores (COLUMN centered data): '
    do i=1,n_measurements_red
       write(1,'(10000f9.3)')(factors_cc(i,j),j=1,n_series_red)
    end do
!
!  Loadings
!
    if(allocated(loadings))deallocate(loadings); allocate(loadings(n_series_red,n_series_red)); loadings=0.d00
    loadings=transpose(eigen_vec_cov)
    loadings=eigen_vec_cov
    write(1,*)
    write(1,*)'Loadings: '
    do i=1,n_series_red
       write(1,'(10000f9.3)')(loadings(i,j),j=1,n_series_red)
    end do
!
!  Scores
!
    if(allocated(scores))deallocate(scores); allocate(scores(n_measurements_red,n_series_red)); scores=0.d00
    scores=matmul(spectra_value_ck,loadings)
    write(1,*)
    write(1,*)'Scores: '
    do i=1,n_measurements_red
       write(1,'(10000f9.3)')(scores(i,j),j=1,n_series_red)
    end do
!
!  Residuals
!
    if(allocated(residuals))deallocate(residuals); allocate(residuals(n_measurements_red,n_series_red));residuals=0.d00
    if(allocated(spectra_odtworzone))deallocate(spectra_odtworzone)
    allocate(spectra_odtworzone(n_measurements_red,n_series_red));spectra_odtworzone=0.d00
    if(allocated(singul))deallocate(singul); allocate(singul(n_measurements_red),vi(n_measurements_red),trv(0:n_measurements_red),erv(n_measurements_red))
    singul=0.d00; erv=0.d00; trv=0.d00
    if(allocated(r8tmp))deallocate(r8tmp); allocate(r8tmp(n_measurements_red,n_series_red));r8tmp=0.d00
    r8tmp=matmul(transpose(scores),scores)
    singul = [(r8tmp(i,i),i=1,n_series_red)]; singul0=sum(singul)
    skw_xauto=sum(spectra_value_ck**2)
    dfij=dfloat(n_series_red*n_measurements_red)
    trv(0)=singul0/dfij
    
    deallocate(r8tmp)
    !write(1,5)n_series_red,(singul(i),i=1,n_series_red)
5   format(/'N= ',i5,' SINGUL matrix (T(transp)*T diagonales) ',/10000f9.3/)
    open(3,file='residuals')
    open(4,file='spectra_odtworzone.dat')
    do i=1,12
        call get_residuals(i)
    end do
    write(1,*)
    write(1,*)' PC     Lambda     TRV      ERV'
    write(1,10)0,singul0,trv(0),0.
10  format(1x,i4,2x,f7.2,2x,f8.3,3x,f8.3)
    do i=1,n_series_red
        write(1,10)i,singul(i),trv(i),erv(i)
    end do
!
!  Factors scores and loadings
!
     if(allocated(factor_loadings))deallocate(factor_loadings)
     allocate(factor_loadings(n_series_red,n_series_red)); factor_loadings=0.d00
     do i=1,n_series_red
        factor_loadings(:,i)=loadings(:,i)*dsqrt(eigvalues_cov(i))
     end do
!    write(1,*)
!    write(1,*)'Factor loadings: '
!    do i=1,n_series_red
!       write(1,'(10000f9.3)')(factor_loadings(i,j),j=1,n_series_red)
!    end do
    if(allocated(factor_scores))deallocate(factor_scores)
    allocate(factor_scores(n_measurements_red,n_series_red)); factor_scores=0.d00
    factor_scores=matmul(selected_spectra_value,factor_loadings)
!    write(1,*)
!    write(1,*)'Factor Scores: '
!    do i=1,n_measurements_red
!       write(1,'(10000f9.3)')(factor_scores(i,j),j=1,n_series_red)
!1    end do
       
   ! factors_rc=matmul(spectra_value_cw,b_matrix)
   ! write(1,*)
   ! write(1,*)'Factor scores (RAW centered data): '
   ! do i=1,n_measurements_red
   !    write(1,'(10000f9.3)')(factors_rc(i,j),j=1,n_series_red)
   ! end do
!
!  Now perform NIPALS procedure
!
    end program file_op
    
    subroutine alloc_eigrs
    use cv
    if(allocated(eigvalues_cov))deallocate(eigvalues_cov)
    allocate(eigvalues_cov(n_series_red));eigvalues_cov=0.d00
    if(allocated(eigen_vec_cov))deallocate(eigen_vec_cov)
    allocate(eigen_vec_cov(n_series_red,n_series_red));eigen_vec_cov=0.d00
    if(allocated(x))deallocate(x)
    n=n_series_red
    allocate(x(n*(n+1)/2+n));x=0.d00
    ier=0
    if(allocated(s))deallocate(s); allocate(s(n_series_red,n_series_red)); s=0.d00
    if(allocated(xtmp_2d))deallocate(xtmp_2d); allocate(xtmp_2d(n_series_red,n_series_red)); xtmp_2d=0.d00
    if(allocated(b_matrix))deallocate(b_matrix); allocate(b_matrix(n_series_red,n_series_red)); b_matrix=0.d00
!
!  factors_nc - Non Centered
!  factors_cc - Column Centered
!  factors_rc - Row Centered
!
    if(allocated(factors_nc))deallocate(factors_nc); allocate(factors_nc(n_measurements_red,n_series_red)); factors_nc=0.d00
    if(allocated(factors_cc))deallocate(factors_cc); allocate(factors_cc(n_measurements_red,n_series_red)); factors_cc=0.d00
    if(allocated(factors_rc))deallocate(factors_rc); allocate(factors_rc(n_measurements_red,n_series_red)); factors_rc=0.d00
    end subroutine alloc_eigrs

   subroutine process_inp
   use cv; use upplow
   character(150) :: a,rows_line,columns_line
   character(1), allocatable :: r(:),c(:)
   logical :: exist
   row_columns_corrections=.false.
   inquire(file='inp.dat',exist=exist)
   if(.not.exist)return
   row_columns_corrections=.true.
   nl=0;rows_line(:)=' '; columns_line(:)=' '
   call get_file('inp.dat',nl,rows_columns,150)
   do i=1,nl
       a=rows_columns(i)
       call low2up(a)
       ipos=index(trim(a),'ROWS:')
       if(ipos/=0)rows_line=a(ipos+5:)
       ipos=index(trim(a),'COLUMNS:')
       if(ipos/=0)columns_line=a(ipos+8:)
   end do
!
   ll=len_trim(rows_line)
   do i=1,ll
     if(rows_line(i:i)==',')rows_line(i:i)=' '
   end do
   if(len_trim(rows_line)==0)goto 1
   call tokenize_string(rows_line)
   ipos=index(string_tokenized(1),'-')
   if(ipos==0)then
     read(string_tokenized(1),*)ii
     if(ii==0)goto 1
   end if
!
   rows=.false.;n_rows=size(rows)
   do i=1,num_tokens
     a(:)=' '
     a=string_tokenized(i)
     ipos=index(a,'-')
     if(ipos==0)then
          read(a,*)ii
          rows(ii)=.true.
       else
         a(ipos:ipos)=' '
         read(a,*)ii1,ii2
         if((ii1>n_rows).or.(ii2>n_rows))then
             write(*,*)'Error in INP1 file - too many rows !!!'
             write(*,*)'ii1,ii2,N_ROWS=',ii1,ii2,n_rows
             stop 'error termination !!!'
         end if
         rows(ii1:ii2)=.true.
     end if
    end do
!
!  Now columns
!
1 continue
   ll=len_trim(columns_line)
   do i=1,ll
     if(columns_line(i:i)==',')columns_line(i:i)=' '
   end do
   if(len_trim(columns_line)==0)return
   call tokenize_string(columns_line)
   ipos=index(string_tokenized(1),'-')
   if(ipos==0)then
     read(string_tokenized(1),*)ii
     if(ii==0)return
   end if
!
   columns=.false.
   do i=1,num_tokens
     a(:)=' '
     a=string_tokenized(i)
     ipos=index(a,'-')
     if(ipos==0)then
          read(a,*)ii
          columns(ii)=.true.
       else
         a(ipos:ipos)=' '
         read(a,*)ii1,ii2
         columns(ii1:ii2)=.true.
     end if
    end do
    end subroutine process_inp
    
    subroutine center
    use cv
    implicit real(8) (a-h,o-z)
!
!   Centrowanie kolumn - centrowanie kolumnowe
!
    if(allocated(aves_row))deallocate(aves_row); allocate(aves_row(n_measurements_red));aves_row=0.d00
    if(allocated(aves_column))deallocate(aves_column); allocate(aves_column(n_series_red));aves_column=0.d00
    if(allocated(spectra_value_cw))deallocate(spectra_value_cw); allocate(spectra_value_cw(n_measurements_red,n_series_red))
    if(allocated(spectra_value_ck))deallocate(spectra_value_ck); allocate(spectra_value_ck(n_measurements_red,n_series_red))
    if(allocated(sdt_dev_col))deallocate(sdt_dev_col); allocate(sdt_dev_col(n_series_red)); sdt_dev_col=0.d00
    if(allocated(sdt_dev_row))deallocate(sdt_dev_row); allocate(sdt_dev_row(n_measurements_red)); sdt_dev_row=0.d00
    spectra_value_ck=0.d00; spectra_value_cw=0.d00; sdt_dev_col=0.d00

    if(center_data)then
       do i=1,n_series_red
           aves_column(i)=average_table(selected_spectra_value(1:n_measurements_red,i))
           sdt_dev_col(i)=std_dev_sample(selected_spectra_value(1:n_measurements_red,i),n_measurements_red)
           !sdt_dev_col(i)=std_dev_pop(selected_spectra_value(1:n_measurements_red,i),n_measurements_red)
           !spectra_value_ck(1:n_measurements_red,i)=(selected_spectra_value(1:n_measurements_red,i)-aves_column(i))/sdt_dev_col(i)
           !spectra_value_ck(1:n_measurements_red,i)=(selected_spectra_value(1:n_measurements_red,i)-aves_column(i))
           spectra_value_ck(1:n_measurements_red,i)=(selected_spectra_value(1:n_measurements_red,i)-aves_column(i))/(sdt_dev_col(i))
       end do
     else
         spectra_value_ck=selected_spectra_value
    end if
!
!  Write reduced spectra data to disk file SELECTED_SPECTRA.DAT
!
    open(1,file='selected_spectra_ck.dat')
    form(:)=' '
    form(1:1)='('
    write(form(2:),'(i3)')n_series_red
    form = trim(form) // '(1x,a9))'
    write(1,trim(form))series_names_red
!
    form(:)=' '
    form(1:4)='(a,'
    write(form(5:),'(i3)')n_series_red
    form = trim(form) // '(1x,f9.4))'
    do i=1,n_measurements_red
        write(1,form)trim(sample_table_red(i)),spectra_value_ck(i,1:n_series_red)
    end do 
    close(1,status='keep')
!
!  Centrowanie wierszowe
!
    sdt_dev_row=0.d00
    if(center_data)then
         do i=1,n_series_red
            aves_row(i)=average_table(selected_spectra_value(i,1:n_series_red))
            sdt_dev_row(i)=std_dev_sample(selected_spectra_value(i,1:n_series_red),n_measurements_red)
            spectra_value_cw(1:n_measurements_red,i)=  &
            (selected_spectra_value(1:n_measurements_red,i)-aves_row(i))/sdt_dev_row(i)
         end do
      else
         spectra_value_cw= selected_spectra_value
    end if
    open(1,file='selected_spectra_cw.dat')
    form(:)=' '
    form(1:1)='('
    write(form(2:),'(i3)')n_series_red
    form = trim(form) // '(1x,a9))'
    write(1,trim(form))series_names_red
!
    form(:)=' '
    form(1:4)='(a,'
    write(form(5:),'(i3)')n_series_red
    form = trim(form) // '(1x,f9.4))'
    do i=1,n_measurements_red
        write(1,form)sample_table_red(i),spectra_value_cw(i,1:n_series_red)
    end do 
    close(1,status='keep')
    end subroutine center
    
    subroutine process_config_file
    use cv
    character(200) :: aa
!
! Remove comments first 
!
    do i=1,nlines
      call low2up(conf_file(i))
      ipos=index(conf_file(i),'!')
      if(ipos==0)cycle
      conf_file(i)(ipos:)=' '
    end do
!
!  Input file name
!
    do i=1,nlines
      aa=conf_file(i)
      ipos=index(aa,'INPUT_NAME')
      if(ipos==0)cycle
      ipos=index(aa,'=')
      ib=0; ie=0
      SINGLE_LINE_BEGIN: do j=ipos+1,len(aa)
        if(aa(j:j)/=' ')then
          ib=j; ie=j
          exit
        end if
      end do SINGLE_LINE_BEGIN
      SINGLE_LINE_END: do j=ib+1,len(aa)
        if(aa(j:j)==' ')then
          ie=j-1
          exit
        end if
      end do SINGLE_LINE_END
      inp_file_name=aa(ib:ie)
    exit
    end do 
!
!  Center data or not
!
    center_data=.true.
    do i=1,nlines
      aa=conf_file(i)
      ipos=index(aa,'CENTER_DATA')
      if(ipos==0)cycle
      ipos=index(aa,'=')
      ib=0; ie=0
      SINGLE_LINE: do j=ipos+1,len(aa)
        if(aa(j:j)=='Y')then
           exit  
         elseif(aa(j:j)=='N')then
           center_data=.false.
           exit
         end if
      end do SINGLE_LINE
    end do
!
!  Diagonalize covariance or correlation matrix
!
    diagonalize_correlation_matrix=.true.
    do i=1,nlines
      aa=conf_file(i)
      ipos=index(aa,'DIAGONALIZE')
      if(ipos==0)cycle
      ipos=index(aa,'=')
      ib=0; ie=0
      KINDD_LINE: do j=ipos+1,len(aa)
        if(aa(j:j)=='R')then
           exit  
         elseif(aa(j:j)=='V')then
           diagonalize_correlation_matrix=.false.
           exit
         end if
      end do KINDD_LINE
    end do
    end subroutine process_config_file
    
    subroutine get_residuals(n)
    use cv
   ! write(3,*)'**********************************************************************************'
    write(3,*); write(4,*)
    write(3,*)'Residuals, Number of components= ',n
    write(4,*)'Spectra odtworzone, Number of components= ',n
    spectra_odtworzone=matmul(scores(:,1:n),transpose(loadings(:,1:n)))
    residuals=spectra_value_ck-spectra_odtworzone
    do i=1,n_measurements_red
       write(3,'(10000f9.4)')(residuals(i,j),j=1,n_series_red)
        write(4,'(10000f9.4)')(spectra_odtworzone(i,j),j=1,n_series_red)
    end do
    do i=1,n_measurements_red
      vi(i)=sum(residuals(i,:)**2)
    end do
    v0=sum(vi)/dfloat(n_measurements_red)
    trv(n)=(singul0-sum(singul(1:n)))/dfij
    erv(n)=1.d00-trv(n)/trv(0)
    write(3,*)
    write(3,1)
1   format(100('*'))
    end subroutine get_residuals