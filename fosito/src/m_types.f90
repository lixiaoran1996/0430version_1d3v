module m_types
  implicit none
  public

  ! real types
  integer, parameter :: sp = kind(0.0)
  integer, parameter :: dp = kind(0.0d0)

  ! character default lengths
  integer, parameter :: tiny_len = 20
  integer, parameter :: name_len = 40
  integer, parameter :: line_len = 200
  
end module m_types
