integer, parameter :: TS_INIT=1
integer, parameter :: TS_INPUT=2
integer, parameter :: TS_INTEG=3
integer, parameter :: TS_OUTPUT=4
integer, parameter :: TS_POST=5
character(len=8), parameter :: names(5) = (/&
"INIT    ", &
"INPUT   ", &
"INTEGRA ", &
"OUTPUT  ", &
"POST    "/)
