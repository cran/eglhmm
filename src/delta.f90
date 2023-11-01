! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine delta(i,j,ijd)
if(i==j) then
    ijd = 1
else
    ijd = 0
endif
end subroutine delta
