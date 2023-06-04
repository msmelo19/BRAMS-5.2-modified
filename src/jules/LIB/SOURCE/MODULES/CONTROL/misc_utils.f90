!#######################################################################
!#######################################################################
!#######################################################################

! module misc_utils
! Contains various miscellaneous procedures - what a useful description.

!###############################################################################

  MODULE misc_utils

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Interfaces.
!-------------------------------------------------------------------------------

  INTERFACE varValue
  MODULE PROCEDURE varValue1D,varValue2D,varValue3D   !  ripe for rationalisation!
  END INTERFACE

  INTERFACE inList
! In future could add other versions - e.g. for character variables.
  MODULE PROCEDURE inListInt
  END INTERFACE inList

!-------------------------------------------------------------------------------

  INTEGER, SAVE :: undefInt = -999 !  missing data value
  REAL, SAVE :: undef = 1.0e20  !  missing data value

  CONTAINS

!#######################################################################
!#######################################################################

! function inListInt
! Internal procedure in module misc_utils.
! Module procedure with generic name inList.
!
! Check that an integer appears in a list.
! Expected use is to check that the value of an option is acceptable.

  FUNCTION inListInt( val,list,name,errMess )  &
                RESULT( found )

  IMPLICIT NONE

! Scalar function result
  LOGICAL ::   &
    found  !  T means input value found in list, F means not in list

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) :: val   !  the value that is to be checked
  CHARACTER(len=*), INTENT(in) :: errMess  !  a message printed on error
  CHARACTER(len=*), INTENT(in) :: name     !  the name of a variable
!                                               (for diagnostic use only)

! Array arguments with intent(in)
  INTEGER, INTENT(in) :: list(:)   !  a list of acceptable values

! Local scalars.
  INTEGER :: i,n

!-------------------------------------------------------------------------------

! Initialise result.
  found = .FALSE.

! Get size of list.
  n = SIZE( list(:) )

! Search list for the input value.
  DO i=1,n
    IF ( val == list(i) ) THEN
      found = .TRUE.
      EXIT
    ENDIF
  ENDDO

! Report an error.
  IF ( .NOT. found ) THEN
    WRITE(*,*)'ERROR: inListInt: value not in list.'
    WRITE(*,*)'name=',TRIM(name),' value=',val
    IF ( n > 0 ) THEN
      WRITE(*,*) 'Size of list=',n
      WRITE(*,*) 'List values=',list(:)
    ELSE
      WRITE(*,*)'List is empty.'
    ENDIF
    WRITE(*,*) TRIM( errMess )
  ENDIF

  END FUNCTION inListInt
!#######################################################################
!#######################################################################
!#######################################################################
  SUBROUTINE checkVars( nvarMax,allowConst,allowSpecial,varUse,foundVar  &
                       ,varDesc,varName,nvar,varFlag,errFound )

! For use in model initialisation, to check that all required input
! variables have been correctly specified (e.g. given a location in a file).

  IMPLICIT NONE

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) :: nvarMax      !  max possible number of variables
  LOGICAL, INTENT(in) :: allowConst   !  switch to allow a constant value
!                                          to be used to set a field
  LOGICAL, INTENT(in) :: allowSpecial !  switch to allow a field to be
!                                          treated as a "special case"
!                                          with extra code required elsewhere

! Array arguments with intent(in)
  INTEGER, INTENT(in) :: varUse(nvarMax) !  flag indicating if a variable
!                        is required for current setup, and how it is used
!              0: not required
!              1: needed. Can be set indirectly (from other variables).
!              2: only needed indirectly to set other variables.
  LOGICAL, INTENT(in) :: foundVar(nvarMax)  !  flag indicating if variable
!              has been found in list.
  CHARACTER(len=*), INTENT(in) :: varDesc(nvarMax)  !  description of each variable
  CHARACTER(len=*), INTENT(in) :: varName(nvarMax)  !  name of each variable

! Scalar arguments with intent(inout)
  INTEGER, INTENT(inout) :: nvar     !  number of variables to be processed

! Array arguments with intent(inout)
  INTEGER, INTENT(inout) :: varFlag(nvarMax) !  flag indicating how a variable is to
!             be set
!             >1: field number (position in file) to be used
!              0: variable not required
!             -1: variable will be set to a constant value
!           <=-2: variable will be set using extra code provided elsewhere (e.g.
!                  using a function of other variables)
!             If a variable is found, but not needed, varFlag is set to zero on output.

! Scalar arguments with intent(out).
  LOGICAL, INTENT(out) :: errFound   !  error flag

! Local scalar variables.
  INTEGER :: ivar  !  loop counter

!----------------------------------------------------------------------
! Initialise.
  errFound = .FALSE.

! Check that all necessary variables are provided.
  DO ivar=1,nvarMax

    IF ( varUse(ivar) > 0 ) THEN

      IF ( foundVar(ivar) ) THEN

!       Variable was found in the list.
!       Determine source of data.

        SELECT CASE ( varFlag(ivar) )

          CASE ( 1: )
!           A location is given. Nothing more to do.
          CASE ( 0 )
!           Variable is in the list, but we don't know how to get it.
            errFound = .TRUE.
            WRITE(*,*)'ERROR: checkVarPos: varFlag=0'
            WRITE(*,*)'Variable=',TRIM(varName(ivar))
            WRITE(*,*)'Don''t know where to get data.'
          CASE ( -1 )
!           A constant value will be used. Check this is allowed.
            IF ( .NOT. allowConst ) THEN
              errFound = .TRUE.
              WRITE(*,*)'ERROR: checkVarPos: varFlag=-1 but allowConst=F'
              WRITE(*,*)'Variable=',TRIM(varName(ivar))
              WRITE(*,*)'Flag indicates that a constant value is to be used,'
              WRITE(*,*)'but switch is set to NOT allow this.'
            ENDIF
          CASE default
!           This covers values <-1. These flags are used to indicate
!           "special cases" specific to a particular section of the model.
!           Check this is allowed.
            IF ( .NOT. allowSpecial ) THEN
              errFound = .TRUE.
              WRITE(*,*)'ERROR: checkVarPos: varFlag<-1 but allowSpecial=F'
              WRITE(*,*)'Variable=',TRIM(varName(ivar))
              WRITE(*,*)'Flag indicates that this variable is to be treated as'
              WRITE(*,*)'a "special case" with specific code, but the switch'
              WRITE(*,*)'has been set to NOT allow this.'
            ENDIF

        END SELECT

      ELSE

!       NOT foundVar
!       Variable is to be used, but was not in list.
!       This is OK only if varFlag has been set to <-1 (indicating that variable will
!       be derived from others).
        IF ( varFlag(ivar) > -2 ) THEN
          errFound = .TRUE.
          WRITE(*,"(/,2a)")'ERROR: checkVars: variable not supplied: ',TRIM(varName(ivar))
          WRITE(*,*) TRIM(varDesc(ivar))
          WRITE(*,*)'This variable is required for the current model configuration.'
        ENDIF

      ENDIF  !  foundVar

    ELSE

!     varUse=0: variable is not required.
      IF ( foundVar(ivar) ) THEN
        WRITE(*,*) 'checkVars: ',TRIM(varName(ivar))  &
            ,' is provided but is not needed. Ignoring.'
        varFlag(ivar) = 0
!       Decrement number of variables to be set.
        nvar = nvar - 1
      ENDIF

    ENDIF   !  varUse
  ENDDO  !  ivar

  END SUBROUTINE checkVars
!#######################################################################
!#######################################################################
!#######################################################################
! function checkVarPos
! Internal procedure in module misc_utils
!
! Check that a list of positions in which variables are to be found
! contains no repeats and no overlaps, and count number of variables
! indicated. Also check that any mandatory variables are present.
! If flag is set, will
! check that positions are increasing (so no need to backspace) and flag
! may also indicate that variables have to be sequential from start.
!
! A useful extension might be to include number of levels of data
! required for each variable - then we could check that there is no
! overlap between data required for different variables.
!
! The assumption throughout is that we do not want to use the same
! data to set value of more than one variable!

  FUNCTION checkVarPos( orderFlag,varPos,message,nfield,geZero,varNzIn )  &
                RESULT( nvar )

  IMPLICIT NONE

  INTEGER ::   &!  scalar function result
    nvar  !  # of variables counted
!             A value <0 indicates an error.

  INTEGER, INTENT(in) ::  &!  IN scalars
    orderFlag   !  flag indicating whether order of variables in
!                  file (as given by varPos) is important
!    =0  order of variables is unimportant and so not checked
!    =1  variable position must increase (so no need to backspace),
!        but varPos<1 ignored
!    =2  variables must be sequential at start of file but varPos<1 ignored

  INTEGER, INTENT(in), OPTIONAL ::  &!  optional in scalars
    nfield    !  number of fields (levels) in file

  INTEGER, INTENT(in) ::  &!  IN arrays
    varPos(:)  !  list of positions of variables
!                   =0 means a variable is not required (and so is ignored)
!                   >0 means position given is location in file
!                   <0 means a variable will be dervied via other means (and so is ignored)

! Optional array arguments with intent(in)
  INTEGER, INTENT(in), OPTIONAL :: varNzIn(:)  !  number of levels (fields)
!                          required for each variable

  LOGICAL, INTENT(in), OPTIONAL ::  &!  optional in scalars
    geZero   !  TRUE means all varPos must be >= 0

  CHARACTER(len=*), INTENT(in) :: message  !  error message

  INTEGER ::   &!  local
    i,j        &!  loop counters
   ,n          &!  size of list
   ,prevVarPos  !  if orderFlag=1 or 2, location of end of previous variable

  LOGICAL ::  &!  scalars
   errFlag    !  TRUE if an error is found

! Local arrays
  INTEGER :: varNz(SIZE(varPos))  !  number of levels (fields)
!                 required for each variable.

!-------------------------------------------------------------------------------
  IF ( orderFlag<0 .OR. orderFlag>2 ) THEN
    WRITE(*,*)'checkVarPos: orderFlag=',orderFlag,' invalid.'
    WRITE(*,*)'Must be 0,1 or 2.'
    WRITE(*,*) TRIM( message )
    STOP
  ENDIF

! If optional argument giving number of levels is not present, set number
! of levels to 1, which means this function will only check that the
! requested varPos values do not contain repeats.
  IF ( PRESENT( varNzIn ) ) THEN
    varNz(:) = varnzIn(:)
  ELSE
    varNz(:) = 1
  ENDIF

  n = SIZE( varPos)
  nvar = 0
  prevVarPos = 0
  errFlag = .FALSE.

  DO i=1,n

!   If flag set, check that all values are >= 0.
    IF ( PRESENT(geZero) ) THEN
      IF ( geZero .AND. varPos(i)<0 ) THEN
        errFlag = .TRUE.
        WRITE(*,*) 'ERROR: checkVarPos: varPos=',varPos(i),' <0'
        WRITE(*,*) TRIM(message)
      ENDIF
    ENDIF

    IF ( varPos(i) > 0 ) THEN

!     Check location is within range available in this file.
      IF ( PRESENT(nfield) ) THEN
        IF ( varPos(i) + varNz(i) - 1 > nfield ) THEN
          errFlag = .TRUE.       
          WRITE(*,*)'varPos=',varPos(i),' out of range.'
          WRITE(*,*)'A variable starting at this location and requiring ',varNz(i)
          WRITE(*,*)'levels will pass the end of file.'
          WRITE(*,*)'Number of fields in file=',nfield
          WRITE(*,*) TRIM(message) 
        ENDIF
      ENDIF

!     Check for overlap (reuse of data).
!     If optional argument varnzIn was not provided, this will just check
!     that the first level of this variable does not match that of any other.
      DO j=i+1,n
        IF ( varPos(j)>0 .AND.  &
             varPos(i)>=varPos(j) .AND. varPos(i)<=varPos(j)+varNz(j)-1 ) THEN
          errFlag = .TRUE.
          WRITE(*,*)'Data required for variable starting at varPos=',varPos(i)
          WRITE(*,*)'overlaps with that used for variable starting at'
            WRITE(*,*)'varPos=',varPos(j)
          IF ( PRESENT(varNzIn) ) WRITE(*,*) 'with ',varNzIn(j),' levels.'
          WRITE(*,*) TRIM(message)
        ENDIF
      ENDDO

!     Check for increasing position.
      IF ( orderFlag==1 ) THEN
        IF ( varPos(i) <= prevVarPos ) THEN
          errFlag = .TRUE.
          WRITE(*,*) 'ERROR: checkVarPos: varPos=',varPos(i),  &
                     ' is before (or at) previous var at varPos=',prevVarPos
          WRITE(*,*) TRIM(message)
        ENDIF
        prevVarPos = varPos(i)  !  update position
      ENDIF

!     Check for sequential positions at start of file.
      IF ( orderFlag == 2 ) THEN
        IF ( varPos(i)/=prevVarPos+1 ) THEN
          WRITE(*,*) 'ERROR: checkVarPos: varPos=',varPos(i)  &
            ,' not sequential to previous var at position=',prevVarPos
          WRITE(*,*)'Or 1st variable is not at varPos=1.'
          WRITE(*,*) TRIM(message)
          errFlag = .TRUE.
        ENDIF
        prevVarPos = varPos(i) + varNz(i) - 1  !  update position
      ENDIF

      nvar = nvar + 1  !  increment count of good variables

    ENDIF   !  varPos>0

  ENDDO

! Indicate error as nvar<0, and give some more help.
  IF ( errFlag ) THEN
    nvar = -1
    WRITE(*,*)'Have varPos=',varPos(:)
    IF ( orderFlag==1 ) THEN
      WRITE(*,*)'orderFlag=1, meaning that locations of variables must be in ascending order.'
    ELSEIF ( orderFlag == 2 ) THEN
      WRITE(*,*)'orderFlag=2, meaning that locations of variables must be in ascending order.'
      WRITE(*,*)'and sequential from start of file (at field=1).'
    ENDIF
  ENDIF

  END FUNCTION checkVarPos
!#######################################################################
!#######################################################################
!#######################################################################
! function repeatVal
! Internal procedure in module misc_utils..
! Checks whether a vector character variable contains repeats.
! Empty values are NOT treated differently (i.e. count as repeats).

  FUNCTION repeatVal( list ) RESULT( resultVal )

  IMPLICIT NONE

  LOGICAL  ::   &!  scalar function result
    resultVal  !   T means there is one or more repeated value in the list
!                  F means there are no repeats

  INTEGER ::  &!  local scalars
    i,j   &!  work
   ,n      !  size

  CHARACTER(len=*), INTENT(in) ::  &!  in arrays
    list(:)    !  input list
!-------------------------------------------------------------------------------

  resultVal = .FALSE.
  n = SIZE( list(:) )

  DO i=1,n-1
    DO j=i+1,n
      IF ( list(i)==list(j) ) THEN
        WRITE(*,*) TRIM(list(i)),' is repeated at locations #',i,' and ',j
        resultVal = .TRUE.
      ENDIF
    ENDDO
  ENDDO

  END FUNCTION repeatVal
!###############################################################################
!###############################################################################
!###############################################################################

! function getWord
! Internal procedure in module misc_utils
! Given a character variable, returns the Nth word (delimited by blanks).
! If the Nth word does not exist, an empty string is returned.
! There's bound to be a better way of doing this....

  FUNCTION getWord(line,n) RESULT( word )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    n     !  the number of the word to be found
    
  INTEGER  ::   &!  local SCALARS
    i,i1,i2  &!  loop counters/work
   ,nfound   &!  number of words found so far
   ,xlen      !  work

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    line    !   input line

  CHARACTER(len=LEN(line)) ::   &!   scalar function result
    word    !  the selected word
!-------------------------------------------------------------------------------
  nfound = 0
  word=''
  i1 = 1
  i2 = 0
  xlen = LEN_TRIM( line )

! Simply loop through the string, looking for a non-blank followed by a blank
! (or a non-blank at end of string).
  DO i=1,xlen
    IF ( line(i:i) /= ' ' ) THEN
      IF ( (i<xlen .AND. line(i+1:i+1)==' ') .OR. i==xlen ) THEN
        nfound = nfound + 1
        i1 = i2 + 1
        i2 = i
        IF ( i==xlen-1 ) i2=xlen
        IF ( nfound == n ) THEN
!         This is the required word.
!         adjustl removes leading blanks if there is more than one blank between words.
          word = ADJUSTL( line(i1:i2) )
          EXIT
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  END FUNCTION getWord
!#######################################################################
!#######################################################################
! subroutine allocate_error
! Internal procedure in module misc_utils
! Prints message warning of error from allocate or deallocate statement.

  SUBROUTINE allocate_error( callType,ierr,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in scalars
    ierr   !   error code from allocate or deallocate

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType  &!  cause of call to this routine - expect to be "alloc" or "dealloc"
   ,errMess    !  message that is printed
!-------------------------------------------------------------------------------

  WRITE(*,*)'WARNING: ',TRIM(errMess)
  IF ( callType == 'alloc' ) THEN
    WRITE(*,*)'WARNING: error from an attempt to use allocate'
  ELSE  !  assume 'dealloc'!
    WRITE(*,*)'WARNING: error from an attempt to use deallocate'
  ENDIF
  WRITE(*,*) 'Error (or sum of errors)=',ierr

  END SUBROUTINE allocate_error
!#######################################################################
!###############################################################################
!###############################################################################

! subroutine getFields

! v4: option to limit number of fields read (switches off error if more found)
! v3: optionally either gets all fields or the field with given number
! v2: deals with "unmatched" quotes.

! Given a character variable, splits it into fields delimited by a given character.
! Anything that is quoted is considered to be a single field, even if it contains the
! delimiter. If the quotes completely enclose a field (i.e. the leading/trailing quote
! is preceded/followed by the delimiter), they are then removed from the field.
! Whitespace is not preserved, except within quotes.
! Because leading and trailing blanks are stripped (leading ones needn't be, but generally
! speeds up), when using " " as the delimiter we can't have empty fields at the start/end
! of a line - they're ignored.
! Results depend partly upon the lengths of the character variables as declared in the calling
! program - if too short fields are truncated.
! There's bound to be some combination of circumstances that I haven't considered and that cause catastrophic failure....
! Note that if the input line has not been initialised (by the calling program), len_trim
! seems to return len(line), and size(field) empty fields are returned.
! Uses assumed-shape arrays, so needs an explicit interface (eg via module).
!-------------------------------------------------------------------------------
! v2: Can deal with "unmatched" quotes.
!     "Matched" quotes are those that appear as the first and last character of a field. Everything inside these, including
!     any occurrence of the delimiting character, is considered to be one field. Other quotes (matched or unmatched) can appear
!     within these matched quotes.
!     May have unnecessary complications that haaven't been removed since earlier versions/have crept in - haven't checked!
! Examples:
!     Input                  delimiter allowRepeat  nfield       fields                                notes
!  23.4, 45.7, 123.4              ,          F         3     #1 23.4  #2 45.7  #3 123.4
!  23.4 | 45.7 | 123.4            |          F         3     #1 23.4  #2 45.7  #3 123.4
!  23.4 || 45.7 || 123.4          |          T         3     #1 23.4  #2 45.7  #3 123.4
! ,23.4, 45.7, 123.4,             ,          F         4     #1       #2 23.4  #3 45.7  #4 123.4    First character is the delimiter.
!    Note that blanks at start of line are ignored, so if the delimiter is ' ', these are ignored
!    at start of line.
! Hello, I'm Big Eck.          blank         F         4     #1 Hello #2 I'm   #3 Big   #4 Eck.
! He said, "I'm Big Eck.".     blank         F         3     #1 He    #2 said, #3 "I'm Big Eck.".
! What,, are you doing?           ,          T         2     #1 What  #2 are you doing?
! Oa, Ob, O'Neill                 ,          F         3     #1 Oa    #2 Ob    #3 O'Neill
! 'He said 'move over'.'       blank         F         1     #1 He said 'move over'.                Note that the enclosing quotes have been stripped.
! 'He said "move over".' O yes. blank        F         3     #1 He said "move over".  #2 O  #3 yes.
! 'He said "O'Phee, move over".' blank       F         1     #1 He said "O'Phee, move over".        Quotes, including unmatched quote, inside quotes.
!
! ***** What DOESN'T work: *******
!    'O dear me.               blank         F         3     #1 'O    #2 dear  #3 me.               Leading quote causes confusion.
!                                             Could only deal with this if scanned ahead and worked out that there was no
!                                             closing quote...but parsing more complicated examples would become difficult! And quite unnecessary for now.

!-------------------------------------------------------------------------------
! v3 Optionally gets either
!    1) all fields
! or 2) the field with given number
! or 3) just counts the fields.
!-------------------------------------------------------------------------------
! v4: If reading a list of fields, optionally stops looking for fields when
!       has all required.
!-------------------------------------------------------------------------------


  SUBROUTINE getFields( clen,allowRepeat,getAll,getNum,delim,lineIn,errFlag  &
                       ,ifield,nfMaxUse,nfield,field,fieldList )

! Use with keywords for optional arguments.

  IMPLICIT NONE

!--------------------------------------------------
! imported scalars with intent(in)

  INTEGER, INTENT(in) ::  &
    clen      !  length of output character variables
!                  NB Code does not work if we try to avoid needing this argument
!                     and instead use len=len(inLine) in declarations. But len(inline) is
!                     apparently OK sometimes! I suspect problem arises if used with fieldList.

  LOGICAL, INTENT(in) ::  &
    allowRepeat  &!  T if repeated adjacent delimiters are to be treated as a single delimiter
   ,getAll       &!  T means get all fields from the input line (subject to nfMax)
!                    F means only get the requested field
   ,getNum        !  T means just count the number of fields in the line, don't return any fields
!                    F means return the requested field(s)

  CHARACTER(len=1), INTENT(in) ::  &
    delim    !  the character that separates the fields
!            !  Cannot be a quotation mark (" or ').

  CHARACTER(len=*), INTENT(in) ::  &
    lineIn    !  the string that is to be split into fields

!----------------------------------------
! imported scalars with intent(out)

  LOGICAL, INTENT(out) ::  &
    errFlag   !  T means an error was raised in this subroutine

!--------------------------------------------------
! imported OPTIONAL scalars wih intent(in)

  INTEGER, INTENT(in), OPTIONAL ::  &
    ifield   !  the number of the requested field (only if getAll=FALSE)

  LOGICAL, INTENT(in), OPTIONAL ::  &
    nfMaxUse    !  Limits the search to at most nfMax fields per line. Used with getAll=TRUE. 
!                  If set, avoids raising an error if find more than nfMax fields on line.
!                    TRUE means look for at most nfMax fields per line.
!                    FALSE means an error is raised if a line contains more fields than we can store.
!----------
! imported OPTIONAL scalars with inent(out)

  INTEGER, INTENT(out), OPTIONAL ::  &
    nfield   !  the number of fields found


  CHARACTER(len=clen), INTENT(out), OPTIONAL ::  &
    field    !  the requested field (only if getAll=FALSE)

!--------------------------------------------------
! imported OPTIONAL arrays with intent(out)

  CHARACTER(len=clen), INTENT(out), OPTIONAL ::  &
    fieldList(:)    !  the requested field (only if getAll=TRUE)

!--------------------------------------------------
! local scalars

  INTEGER  ::   &
    i,i1,i2,i3  &!  loop counters/work
   ,ipos        &!  work (position in current field)
   ,nf          &!  counter of number of fields found
   ,nfMax       &!  maximum possible number of fields
   ,xlen         !  work

  LOGICAL ::  &
    doubleQuoted   &!  T when current character is inside "significant" double quotes
   ,singleQuoted    !  T when current character is inside "significant" single quotes

  CHARACTER(len=LEN(lineIn)) ::  &
    line       &!  working copy of lineIn
   ,newField    !  a field from line

!-------------------------------------------------------------------------------

! Initialise arguments with intent(out).
  errFlag = .FALSE.
  IF ( PRESENT(field) ) field=''
  IF ( PRESENT(nfield) ) nfield = 0
  IF ( PRESENT(fieldList) ) fieldList(:) = ''

!-------------------------------------------------------------------------------
! Check the arguments, including that any required optional arguments are present.

  IF ( getAll ) THEN
!   All fields will be returned (up to max). getNum must be FALSE.
    IF ( getNum ) THEN
      WRITE(*,*)'ERROR: getFields: conflicting arguments.'
      WRITE(*,*)'Only one of getAll and getNum can be TRUE.'
      errFlag = .TRUE.
    ENDIF
    IF ( .NOT.PRESENT(nfield) .OR. .NOT.PRESENT(fieldList) .OR. .NOT.PRESENT(nfMaxUse) ) THEN
      WRITE(*,*)'ERROR: getFields: optional arguments missing.'
      WRITE(*,*)'Need nfield, fieldList and nfMaxUse.'
      errFlag = .TRUE.
    ENDIF
  ELSE
!   NOT getAll
!   Either one field is returned, or the fields are counted.
    IF ( getNum ) THEN
      IF ( .NOT.PRESENT(nfield) ) THEN
        WRITE(*,*)'ERROR: getFields: optional argument missing.'
        WRITE(*,*)'Need nfield.'
        errFlag = .TRUE.
      ENDIF
    ELSE
      IF ( .NOT.PRESENT(ifield) .OR. .NOT.PRESENT(field) ) THEN
        WRITE(*,*)'ERROR: getFields: optional argument(s) missing.'
        WRITE(*,*)'Need ifield and field.'
        errFlag = .TRUE.
      ENDIF
      IF ( ifield < 1 ) THEN
        WRITE(*,*)'WARNING: getFields: ifield<1 - does not exist.'
        errFlag = .TRUE.
      ENDIF
    ENDIF  !  getNum
  ENDIF  !  getAll

! Abandon, if an error is indicated.
  IF ( errFlag ) RETURN

!-------------------------------------------------------------------------------
! Initialise local variables.
  doubleQuoted = .FALSE.
  singleQuoted = .FALSE.
  i1 = 1
  i2 = 0
  ipos = 0
  nf = 0
  line = lineIn
  IF ( getAll ) nfMax =  SIZE(fieldList)

! If input string is empty, there's nothing more to do.
  IF ( LEN_TRIM(line) == 0 ) RETURN

! Refuse to use quote as the delimiter.
  IF ( delim=='"' .OR. delim=="'" ) THEN
    WRITE(*,*) 'ERROR: getFields: cannot use quotation mark as delimiter.'
    errFlag = .TRUE.
    RETURN
  ENDIF
!-------------------------------------------------------------------------------

! Remove leading blanks from line.  
  line = ADJUSTL( line )
  xlen = LEN_TRIM( line )

! Here, "delimiter" means the delimiting character.
! Loop through the string, looking for occurrences of the delimiter.
  DO i=1,xlen

    ipos = ipos + 1

!   Check for "significant" quotes - these must be the first or last character in a field,
!   that is, immediately preceded/followed by the delimiter.

    IF ( line(i:i)=="'" .AND. .NOT.doubleQuoted ) THEN
      IF ( .NOT.singleQuoted ) THEN
!       Check if this is the start of a field.
        IF ( i==1 .OR. ((ipos==1 .AND. line(i-1:i-1)==delim) ) ) singleQuoted=.TRUE.
      ELSE
!       Check if this is the end of a field.
        IF ( i==xlen .OR. ( i<xlen .AND. line(i+1:i+1)==delim ) ) singleQuoted=.FALSE.
      ENDIF
    ENDIF

    IF ( line(i:i)=='"' .AND. .NOT.singleQuoted ) THEN
      IF ( .NOT.doubleQuoted ) THEN
!       Check if this is the start of a field.
        IF ( (ipos==1 .AND. line(i-1:i-1)==delim) .OR. i==1 ) doubleQuoted=.TRUE.
      ELSE
!       Check if this is the end of a field.
        IF ( (i<xlen .AND. line(i+1:i+1)==delim) .OR. i==xlen ) doubleQuoted=.FALSE.
      ENDIF
    ENDIF

    IF ( line(i:i)==delim .OR. i==xlen ) THEN

!     We have a delimiter, and/or the end of the line (which is effectively a delimiter).

      IF ( i == xlen ) THEN
!       Override quotes - we want to end the field.
        doubleQuoted = .FALSE.
        singleQuoted = .FALSE.
      ENDIF

!     Deal with repeated delimiters, if these are allowed.

      IF ( i<xlen .AND. allowRepeat .AND. i>1 .AND. line(i-1:i-1)==delim ) THEN
        i2 = i  !  this will form starting delimiter for next field
        CYCLE
      ENDIF
      
!     If we are not in quotes, this is the end of a field.
      IF ( .NOT.doubleQuoted .AND. .NOT.singleQuoted ) THEN
        nf = nf + 1
        IF ( getAll .OR. getNum ) nfield = nf
        IF ( getAll ) THEN
          IF ( nfield > nfMax .AND. .NOT.nfMaxUse ) THEN
            WRITE(*,*) 'WARNING: getFields: max number of fields has been exceeded'
            WRITE(*,*)'nfMax=',nfMax
            nfield = nfMax
            errFlag = .TRUE.
            RETURN
          ENDIF
        ENDIF
        i1 = i2 + 1   !   location of start of field
        i2 = i        !   delimiter at end of field
        i3 = i2 - 1   !   end of field
        IF ( i==xlen .AND. line(i:i)/=delim ) i3=xlen
        IF ( i3-i1+1 > clen ) WRITE(*,*)'WARNING: getFields: field too long for variable.'
        newField = ADJUSTL( line(i1:i3) )

!       If first and last characters are quotes, remove them.
        i3 = LEN_TRIM( newField )
        IF ( ( newField(1:1)=="'" .AND. newField(i3:i3)=="'" ) .OR. &
             (newField(1:1)=='"' .AND. newField(i3:i3)=='"') )  &
                newField = newField(2:i3-1)

!       Save this field, if appropriate.
        IF ( .NOT. getNum ) THEN
          IF ( getAll ) THEN
            fieldList(nfield) = newField
          ELSEIF ( nf == ifield ) THEN
!           This is the required field.
            field = newField
            EXIT
          ENDIF
        ENDIF

!       If we have found all fields required, return.
        IF ( getAll .AND. nfMaxUse .AND. nfield==nfMax ) RETURN

        ipos = 0   !  ready for next field

      ENDIF  !  quoted
    ENDIF  !   delimiter

  ENDDO

  END SUBROUTINE getFields

!################################################################################
!################################################################################
! subroutine init_count
! Internal procedure in module misc_utils.
! Initialise a counter to be in phase with days.
! e.g. model timestep=1hr, routing timestep=2hr, run starts 01H, initialise counter
! so that routing is done at 02H, 04H,....
! This is based on init_out_count_sub1 - look there if you need more code for other cases.


  SUBROUTINE init_count( period,startDate,startTime,step )

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_360

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    dateToBits

  USE timeConst, ONLY : &
    iSecInDay

  IMPLICIT NONE
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) ::  &!  in scalars
    startDate   &!  date (yyyymmdd) at start of this "section" of run
   ,startTime   &!  time of day (s) at start of this "section" of run
   ,period           !  period (timestep) for process (number of model timesteps)

  INTEGER, INTENT(out) ::  &!  out scalars
    step          !   counter

  INTEGER ::  &!  local SCALARS
    dateDay,dateMonth,dateYear    &!  work
   ,periodDay                      !  period, expressed as number of days

!----------------------------------------------------------------------------

! Set counter so that step=period at the end of the first interval (when process is first activated).
! e.g. 3 hourly process with start time 04H, first call to process will be at 06H
! e.g. 2 day process with start day 1, first call to process will be at OH on day 3 (i.e. end of day 2)
! This is done so that calls to process are in phase with days - this is generally not essential, but
! it's tidier (and means process calls are in phase with output - unless output code is changed!).
! eg 1 hour process timestep, and the run starts at say 15H, we ensure the calls to process are
! in phase with days by initialising the counter so that the calls to process are at 2H, 4H, 6H etc.
! This may mean that the first call to the process occurs before the expected time has elapsed
! (e.g. run starts 02H, 1hr timestep, 2hr process timestep, first call to process will be after only
! 1hr, at 02H), which should be dealt with elsewhere (if important).

  IF ( period > 0 ) THEN

!   First, work out when the process will first be called.

    IF ( period*NINT(timeStep) <= iSecInDay ) THEN

!      Period is less than or equal to one day.
!      One day is known to be a multiple of the period.
!      The first process timestep is complete at the first time in the
!      current day (and after current time) when the time of day is a multiple of the process period
!      (or 00H on next day if period=1day).
!      Add 1 below so that we call the process on the timestep that ENDS at the required time.
       step = MOD( startTime, period*NINT(timeStep) ) / NINT(timeStep) + 1
       IF ( period == 1 ) step = 1
       WRITE(*,*)'init_count startTime=',startTime,' period=',period,' step=',step

    ELSE

!     Period is more than one day (and is previously tested to ensure it is <=30days).
!     Period is known to be a multiple of one day.
!     Generally just call this process starting at 0H on next day
!     [note that this may not match up with output times - see output code],
!     but if using 360-day calendar and period is a factor of 30, keep in phase with months.

      periodDay = period * NINT(timeStep) / iSecInDay
      IF ( l_360 .AND. MOD(30,periodDay)==0 ) THEN
!       e.g. period=2days, call after 2,4,6,... days of month.
!       Get day of month from startDate.
        CALL dateToBits( startDate,dateDay,dateMonth,dateYear,l_360,'init_count')
!       Get mod( timesteps into month at start date/time, period ).
!       Add 1 so that we call the process on the timestep that ENDS at the required time.
        step = MOD( ((dateDay-1)*iSecInDay + startTime)/NINT(timestep), period ) + 1
        WRITE(*,*)'Period is special (l_360 etc), step=',step
      ELSE
!       Work out how many "main" timesteps until 00H on next day. Start time will fall on a timestep.
!       Subtract from period.
!       Add 1 so that we call the process on the timestep that ENDS at the required time.
        step = period - ( iSecInDay - startTime ) / NINT(timestep) + 1
        WRITE(*,*)'Period is days but not special, step=',step
      ENDIF


    ENDIF !  period

!   Subtract one, so that the first increment will recover the value set above.
    step = step - 1

  ELSE

!   period < 0
    WRITE(*,*)'ERROR: init_count: no code for period<=0'
    STOP

  ENDIF

  END SUBROUTINE init_count
  

!################################################################################
!################################################################################
! function replaceTemplate
! Internal procedure in module misc_utils.
! Given a string, generate an output string by replacing any template elements.
! Expected use is to get names of data files.

  FUNCTION replaceTemplate( inString,doTime,errMess,date,time,replString )  &
                    RESULT( outString )

! This function will replace any occurences of search strings (templString)
! with bits of one or more input arguments.
! If doTime=TRUE, date and time are used to provide replacement strings.
! If doTime=FALSE, replString is used as the replacement string.
! If no search strings are found, the result is the input string.
!-------------------------------------------------------------------------------

  USE inout, ONLY :  &
!  imported scalar parameters
     ntemplString,templString  &
!  imported array parameters
    ,metOChar

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_360

  USE time_mod, ONLY :  &
!  imported procedures
     dateToBits,monthName3,secToSMH

  IMPLICIT NONE

!--------------------------------------------------
! Scalar arguments with intent(in)

  CHARACTER(len=*), INTENT(in) ::  &
     inString  &!  an input string (e.g. file name with template strings)
    ,errMess    !  a message to be printed on error

  LOGICAL, INTENT(in) ::  &
    doTime   !  flag indicating what kind of template strings are to be replaced
!                 T means replace time template strings
!                 F means replace other template strings

!--------------------------------------------------
! Optional scalar arguments with intent(in)

  INTEGER, INTENT(in), OPTIONAL ::  &
    date    &!  a date (yyymmdd) to be used to replace time template strings
   ,time     !  a time of day (s) to be used to replace time template strings

  CHARACTER(len=*), INTENT(in), OPTIONAL ::  &
     replString      !  a replacement string to be included in output if doTime=F

!--------------------------------------------------
! Function result
  CHARACTER(len=LEN(inString)) :: outString
! Note that we cannot use an assumed length variable for this result, which means that if
! the required result is longer than len(inString), we end with an error.

!--------------------------------------------------
! Local scalars

  INTEGER ::  &
    dateDay,dateMonth,dateYear  &!  bits of date
   ,h           &!  hour of day
   ,i,i1,i2,ipos   &!  work
   ,j           &!  work
   ,m           &!  minutes of hour
   ,outlen      &!  trimmed length of output string
   ,outlenMax   &!  length of output string
   ,rlen        &!  length of a replacement string
   ,s           &!  seconds
   ,tlen         !  trimmed length of a template string

  LOGICAL ::  &
    found    !  work

  CHARACTER(len=LEN(templString)) ::  &!  local scalars
    tString    !  used to hold a template string

  CHARACTER(len=100) ::  &
    workChar   !  used to hold replacement for template string.
!                   Must be long enough for longest possible replacement.

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Check for optional arguments.
  IF ( doTime ) THEN
    IF ( .NOT. ( PRESENT(date) .AND. PRESENT(time) ) ) THEN
      WRITE(*,*)'ERROR: replaceTemplate: optional argument(s) missing'
      WRITE(*,*)'doTime=T requires date and time to be passed'
      WRITE(*,*) TRIM(errMess)
      STOP
    ENDIF
  ELSE
    IF ( .NOT. ( PRESENT(replString) ) ) THEN
      WRITE(*,*)'ERROR: replaceTemplate: optional argument missing'
      WRITE(*,*)'doTime=F requires replacement string to be passed.'
      WRITE(*,*) TRIM(errMess)
      STOP
    ENDIF
  ENDIF

! Initialise result.
  outString = inString

! Initialise local variables.
  outlen = LEN_TRIM( outString )
  outlenMax = LEN( outString )

! Split date.
  IF ( doTime ) CALL dateToBits( date,dateDay,dateMonth,dateYear,l_360,errMess )

! Go through the string, searching for template strings.
! Note that these are all assumed to be the same length.
  i1 = 1
  DO
    IF ( i1 > outlenMax ) EXIT

!   Search for a %.
    ipos = INDEX( outString(i1:outlen), '%' )
    IF ( ipos == 0 ) EXIT
    ipos = ipos + i1 - 1   !  location in outString
    i2 = ipos + LEN(templString) - 1
    IF ( i2 > outLenMax ) THEN
      WRITE(*,*)'ERROR: replaceTemplate: character variable not long enough.'
      WRITE(*,*)'Length of outString=',LEN(outString)
      WRITE(*,*)'Currently need at least ',i2
      STOP
    ENDIF
    tString = outString(ipos:i2)

!   Check we have a recognised template string.
    found = .FALSE.
    DO i=1,ntemplString
      IF ( tString == templString(i) ) found = .TRUE.
    ENDDO

    IF ( .NOT. found ) THEN
      WRITE(*,*)'ERROR: replaceTemplate: do not recognise template string: ',tString
      WRITE(*,*)'inString=',TRIM(inString)
      WRITE(*,"(a,(15(tr1,a)))")'Valid template strings:',(TRIM(templString(j)),j=1,ntemplString)
      WRITE(*,*)TRIM( errMess )
      STOP
    ENDIF

!   Get trimmed length of template string.
    tlen = LEN_TRIM( tstring )

!   If this template string is not relevant for this call, it is not replaced,
!   but we have to update position for next search.
    SELECT CASE ( tstring )
      CASE ( '%tc','%y2','%y4','%yc','%m1','%m2','%mc','%mm'  &
            ,'%d1','%d2','%dc','%h1','%h2','%hc','%n2' )
!       Time strings.
        IF ( .NOT. doTime ) THEN
          i1 = ipos + tlen
          CYCLE
        ENDIF
      CASE ( '%vv' )
!       Not a time string..
        IF ( doTime ) THEN
          i1 = ipos + tlen
          CYCLE
        ENDIF
      CASE default
        WRITE(*,*)'ERROR: replaceTemplate: no code for template string=',tString
        STOP
    END SELECT

!   Replace template string.
    SELECT CASE ( tstring )
      CASE ( '%tc' )
!       Decades since 1800.
        rlen = 1
!       This is cyclic (as apparently used in UM), so avoid very long runs!
        workChar = metOChar( MOD( (dateYear-1800)/10, SIZE(metOChar) ) )
      CASE ( '%yc' )
!       Year in decade.
        rlen = 1
        workChar = metOChar( MOD( dateYear, 10 ) )
      CASE ( '%y4' )
        rlen = 4
        WRITE(workChar,"(i4.4)") dateYear
      CASE ( '%y2' )
        rlen = 2
        WRITE(workChar,"(i2.2)") MOD(dateYear,100)
      CASE ( '%mm' )
        rlen = 1
        workChar = metOChar( dateMonth )
      CASE ( '%m2' )
        rlen = 2
        WRITE(workChar,"(i2.2)") dateMonth
      CASE ( '%m1' )
        IF ( dateMonth < 10 ) THEN
          rlen = 1
          WRITE(workChar,"(i1.1)") dateMonth
        ELSE
          rlen = 2
          WRITE(workChar,"(i2.2)") dateMonth
        ENDIF
      CASE ( '%mc' )
        rlen = 3
        workChar(1:3) = monthName3( dateMonth )
      CASE ( '%dc' )
        rlen = 1
        workChar = metOChar( dateDay )
      CASE ( '%d2' )
        rlen = 2
        WRITE(workChar,"(i2.2)") dateDay
      CASE ( '%d1' )
        IF ( dateDay < 10 ) THEN
          rlen = 1
          WRITE(workChar,"(i1.1)") dateDay
        ELSE
          rlen = 2
          WRITE(workChar,"(i2.2)") dateDay
        ENDIF
      CASE ( '%hc' )
        rlen = 1
        CALL secToSMH( time,s,m,h,errMess )
        workChar = metOChar( h + 1 )
      CASE ( '%h2' )
        rlen = 2
        CALL secToSMH( time,s,m,h,errMess )
        WRITE(workChar,"(i2.2)") h
      CASE ( '%h1' )
        rlen = 2
        CALL secToSMH( time,s,m,h,errMess )
        IF ( h < 10 ) THEN
          rlen = 1
          WRITE(workChar,"(i1.1)") h
        ELSE
          rlen = 2
          WRITE(workChar,"(i2.2)") h
        ENDIF
      CASE ( '%n2' )
        rlen = 2
        CALL secToSMH( time,s,m,h,errMess )
        WRITE(workChar,"(i2.2)") m
      CASE ( '%vv' )
        rlen = LEN_TRIM(replString)
        IF ( rlen > LEN(workChar) ) THEN
          WRITE(*,*)'ERROR: replaceTemplate: work space not long enough.'
          WRITE(*,*)'rlen=',rlen,' len(workChar)=',LEN(workChar)
          WRITE(*,*)'Increase length of workChar.'
          WRITE(*,*)(TRIM(errMess))
          STOP
        ENDIF
        workChar = TRIM(replString)
      CASE default
        WRITE(*,*)'ERROR: replaceTemplate: no code for template string=',tString
        STOP
    END SELECT

!   Check there is enough space for replacement.
    IF ( outlen+rlen-tlen > outlenMax ) THEN
      WRITE(*,*)'ERROR: replaceTemplate: character variable not long enough.'
      WRITE(*,*)'Length of outString=',LEN(outString)
      WRITE(*,*)'Currently need at least ',outLen+rlen-tlen
      WRITE(*,*)'outString=',TRIM(outString)
      STOP
    ENDIF

!   Replace template string.
    outString = outString(1:ipos-1) // workChar(1:rlen) // outString(ipos+tlen:)

!   Update length of output string.
    outlen = LEN_TRIM( outString )

!   Set index for next loop.
!   We go to position just after the just-replaced template string (i.e.
!   ipos+len(templString)), then adjust for any difference between the lengths
!   of the template and replacement strings (i.e. len_templString)-repl).
!   This gives ipos+len(templString)-(len_templString)-repl), which simplifies
!   to ipos+repl.
    i1 = ipos + rlen

  ENDDO   !  loop through input string

  END FUNCTION replaceTemplate
!###############################################################################
!###############################################################################

! function rm_path
! Internal procedure in module misc_utils.
! Given a filename (incl path), remove all directory from path.

  FUNCTION rm_path( inPath ) &
                    RESULT( fileName )

!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Scalar arguments with intent(in)
  CHARACTER(len=*), INTENT(in) ::  &
     inPath   !  an input string (expected to be file name)

! Function result
  CHARACTER(len=LEN(inPath)) :: fileName

! Local scalars
  INTEGER :: ipos

!--------------------------------------------------

! Initialise result.
  fileName = ''

! Look for last "/" in path.
  ipos = INDEX( inPath, '/', .TRUE. )
  IF ( ipos > 1 ) THEN
    fileName = inPath(ipos+1:)
  ELSE
    fileName = inPath
  ENDIF

  END FUNCTION rm_path
!###############################################################################
!###############################################################################

! subroutine find_file
! Internal procedure in misc_utils.
! Given lists of file names and start times, work out what file the given time
! lies in. Also work out what time level (timestep number) within the file the
! requested time represents (return time at or before requested time).

  SUBROUTINE find_file( time,date,climatol,endTime,next,dataPer,itimeStep  &
                       ,ifile,t  &
                       ,templateT,filePer,fileTime,fileDate  &
                       ,fileName,fileNameTemplate  &
                       ,templateDate,templateTime,templateUnits,errMess )

  USE inout, ONLY :  &!
!   imported scalar parameters
     periodMon

  USE switches, ONLY :  &!
!   imported scalars with intent(in)
     l_360

  USE time_mod, ONLY :  &!
!  imported procedures
     dateToBits,secToSMH,timeDate,timeDate_cmp,timeDate_diff

  USE timeConst, ONLY : &!
     iSecInDay, iSecInHour, iSecInMin

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  In scalars
    date     &!  the date (yyyymmdd) for which data are required
   ,dataPer  &!  the period of the data in the files
!                  >0 indicates period as a number of timesteps
!                 <=0 indicate "special" cases (e.g. monthly).
   ,filePer  &!  the period of the files
!                  >0 indicates period as a number of timesteps
!                 <=0 indicate "special" cases (e.g. monthly).
   ,templateDate  &!  the date (yyyymmdd) associated with template files
   ,templateTime  &!  the time of day (s) associated with template files
   ,time     &!  the time of day (s) for which data are required
   ,itimeStep  !  timestep length (s) for use with dataPer and filePer

  INTEGER, INTENT(inout) ::  &!  inout arrays
    fileDate(:)   &!  start date (yyyymmdd) of each file,in chronological order
   ,fileTime(:)    !  start time (seconds from 00H) of each file,in chronological order

  INTEGER, INTENT(inout) ::  &!  Inout scalars
    ifile   &!  index (number) of file
!                 in: the index of the file that was last used (used if next=TRUE)
!                 out: the index of the file corresponding to required time and date
   ,t         !  time level within file corresponding to required time and date
!                  On output, if ifile=1 but t<1, this indicates that the time was before the
!                  start of the first file.

  INTEGER ::  &!  Local scalars
    dateDay,dateMonth,dateYear   &!  bits of date
   ,fileDateName  &!  date as used in file naming
   ,fileTimeName  &!  date as used in file naming
   ,h        &!  a number of hours
   ,i        &!  loop counter
   ,inc      &!  increment
   ,iyr      &!  work
   ,lookDate &!  date within year (mmdd)
   ,m        &!  a number of minutes
   ,nfile    &!  length of lists (number of files)
   ,ndy,nhr,nmin,nmon,nsec,nyr   &!  work
   ,replDate  &!  date (yyyymmdd) used to replace time template strings
   ,replTime  &!  time of day (s) used to replace time template strings
   ,s         &!  a number of seconds
   ,templateDay,templateH,templateM,templateMonth,templateS,templateYear  &!  work
   ,tmpDate1,tmpDate2,tmpTime    !  work

  LOGICAL, INTENT(in) ::  &!  in scalars
    climatol  &!  TRUE means files are to be treated as climatological (independent of year)
   ,endTime   &!  TRUE means file names use end time (time of last data)
!                 FALSE means names use time of first data.
!                 Only used with time-templating.
   ,next      &!  TRUE means assume that the required time is the next in the current file,
!                   or is in the next file  
               !  FALSE means search all files for required time
   ,templateT  !  TRUE means file names follow a template, with the time of each file
!                   shown in its name

  CHARACTER(len=*), INTENT(inout) ::  &!  inout arrays
    fileName(:)   !  list of file names

  LOGICAL ::  &!  local scalars
    findT         &!  TRUE means calculate time level within file.
   ,found          !  work

  CHARACTER(len=*), INTENT(in) ::  &!  In scalars
    errMess       &!  message printed on error
   ,fileNameTemplate   &!  template file name
   ,templateUnits       !  units used to describe template

  CHARACTER(len=3) ::  &!  local scalars
    units    !  units for increment

!-------------------------------------------------------------------------------
! Note: At present, JULES always has nfile>1 (so as to allow easy comparison
! with time of next file, ifile+1). nfile=1 is coded for at some points, but not
! all, which is OK for now.
!-------------------------------------------------------------------------------

! Check that dataPer is one that is coded for.
  IF ( dataPer<0 .AND. dataPer/=periodMon ) THEN
    WRITE(*,*)'ERROR: find_file: dataPer<=0 must be periodMon(=',periodMon,')'
    STOP
  ENDIF

! Get size of lists (assume all have same length).
  nfile = SIZE( fileTime )
  IF ( nfile == 0 ) THEN
    WRITE(*,*)'ERROR: find_file: list of files is empty.'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

! Initialise flag.
  findT = .FALSE.

!-------------------------------------------------------------------------------
! Find the file that contains data for the requested time.
!-------------------------------------------------------------------------------

! Note that although the requested time is known to be "consistent" with data
! period (i.e. it is a time that would match exactly with a data time), it may
! in fact be outside the times with available data.

  IF ( next ) THEN

!-------------------------------------------------------------------------------
!   Use the next available time of data.
!-------------------------------------------------------------------------------
    t = t + 1

!   Check if the next time is stored in same file as was last used.
!   Nothing to do if there is only one file (nfile=1), or if we are already using last file.
!   Note that, without knowing end time of file, cannot test if we are beyond end of last file.
!   At present we don't need a separate test on templateT, because nfile=2 and ifile=1, but
!   leaving in for clarity.
    IF ( templateT .OR. (nfile>1 .AND. ifile<nfile) ) THEN

      IF ( timeDate_cmp(time,date,'>=',fileTime(ifile+1),fileDate(ifile+1)  &
          ,'find_file') ) THEN

!-------------------------------------------------------------------------------
!       Data are stored in next file.
!-------------------------------------------------------------------------------
        IF ( .NOT. templateT ) THEN

          ifile = ifile + 1
          t = 1   ! consistent with assumption that next data are used, next data are first in next file
!         Deal with climatological files.
!         Date of last file has been set to be 1 year after date of first file.
          IF ( climatol .AND. ifile==nfile ) THEN
!           Reuse first file.
            ifile = 1
!           Reset dates of all files to be 1 year later.
!           Note that this will fail if we try to access 29th Feb in a non leap year.
            DO i=1,nfile
              CALL  timeDate( fileTime(i),fileDate(i),12,'mon',l_360  &
                          ,tmpTime,tmpDate1,'find_file' )
              fileDate(i) = tmpDate1
            ENDDO
!           We may not be at the start of this file, so set flag so that we recalculate
!           time level in file. Note that we have already assumed that the given time is
!           a time with data.
!           (At present, this does not mean we can use data with period that is not a factor of
!            1 year - e.g. 10days. e.g. data starting 01Jan, 10day period, data are read at days
!            360 and 370. At day 370 we detect that we need to go back, but because find_file
!            does not set dataStep, we do not get correct data.)
            findT = .TRUE.
          ENDIF

        ELSE

!         templateT
!         Update file details.
!         Move 'next' file into current file.
          fileTime(1) = fileTime(2)
          fileDate(1) = fileDate(2)
          fileName(1) = fileName(2)
!         Reset ifile and t.
          ifile = 1
          t = 1   ! consistent with assumption that next data are used, next data are first in next file
!         Get details of new next file, by advancing time by one file period.
!         Set units and increment for file period.
          CALL get_filePer_vars( filePer,itimeStep,inc,units )
!         Get date and time.
          CALL timeDate( fileTime(1),fileDate(1),inc,units,l_360,fileTime(2),fileDate(2),'find_file' )

!         Get date and time as used in name.
          IF ( endTime ) THEN
!           Advance by another file period.
            CALL timeDate( fileTime(2),fileDate(2),inc,units,l_360  &
                          ,fileTimeName,fileDateName,'find_file' )
          ELSE
            fileTimeName = fileTime(2)
            fileDateName = fileDate(2)
          ENDIF

!         Get name of file by substituting in template. 2nd argument (doTime) is TRUE to
!         indicate that we are replacing time template strings only.
          fileName(2) = replaceTemplate( fileNameTemplate,.TRUE.,'find_file'  &
                                        ,fileDateName,fileTimeName )
        ENDIF  !  templateT

      ENDIF  !  next file
    ENDIF   !  templateT or nfile etc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  ELSE   !  NOT next
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!   Start from scratch, and work out what file is needed.

!   Indicate that we will calculate time level within file.
    findT = .TRUE.

!-------------------------------------------------------------------------------
    IF ( templateT ) THEN

      ifile = 1

!     We know the time/date for which we want data. Now work out what time/date
!     to use in file name.
!     For annual, monthly and daily template files (as described by templateUnits),
!     we have insisted that there is one file per year/month/day. For shorter
!     periods (templateUnits)(call the period Tau), files can hold data for
!     n*Tau, n>=1.

!     Set units and increment for file period.
      CALL get_filePer_vars( filePer,itimeStep,inc,units )

!     Get components of time and date from template time and requested time.
      CALL dateToBits( templateDate,templateDay,templateMonth,templateYear,l_360,errMess )
      CALL secToSMH( templateTime,templateS,templateM,templateH,errMess )
      CALL dateToBits( date,dateDay,dateMonth,dateYear,l_360,errMess )
      CALL secToSMH( time,s,m,h,errMess )

!     Assume that we will use template time.
      replTime = templateTime

      SELECT CASE ( templateUnits )
        CASE ( 'yr' )
!         Replace year in template date.
          replDate = dateYear*10000 + templateMonth*100 + templateDay
        CASE ( 'mo' )
!         Replace month and year in template date.
          replDate = dateYear*10000 + dateMonth*100 + templateDay
        CASE ( 'dy' )
!         Replace day, month and year in template date.
          replDate = dateYear*10000 + dateMonth*100 + templateDay
        CASE ( 'hr' )
!         Go to n*filePer H on current date.
          replDate = date
!         Calculate time of day (s) at the required file time.
          replTime = FLOOR( REAL(s)/(REAL(filePer)*itimeStep) ) * filePer * itimeStep
        CASE ( 'mn' )
!         Go to n*filePer mins on current date.
          replDate = date
!         Calculate time of day (s) at the required file time.
          replTime = FLOOR( REAL(s)/(REAL(filePer)*itimeStep) ) * filePer * itimeStep
        CASE default
          WRITE(*,*)'ERROR: find_file: ',TRIM(errMess)
          WRITE(*,*)'No code for templateUnits=',templateUnits
          STOP
      END SELECT

!     Load file time and date, and get file name from template.
      fileTime(1) = replTime
      fileDate(1) = replDate
      fileDateName = fileDate(1)
      fileTimename = fileTime(1)
!     If files are named using end time, account for this.
      IF ( endTime ) CALL timeDate( fileTime(1),fileDate(1),inc,units,l_360  &
                                   ,fileTimeName,fileDateName,'find_file' )
!     Get name of file by substituting in template. 2nd argument (doTime) is TRUE to
!     indicate that we are replacing time template strings only.
      fileName(1) = replaceTemplate( fileNameTemplate,.TRUE.,'find_file'  &
                                    ,fileDateName,fileTimeName )

!     Get details of next file, by advancing time by one file period.
      CALL timeDate( fileTime(1),fileDate(1),inc,units,l_360  &
                    ,fileTime(2),fileDate(2),'find_file' )
      fileDateName = fileDate(2)
      fileTimeName = fileTime(2)
!     If files are named using end time, account for this.
      IF ( endTime ) CALL timeDate( fileTime(2),fileDate(2),inc,units,l_360  &
                                   ,fileTimeName,fileDateName,'find_file' )
!     Get name of file by substituting in template. 2nd argument (doTime) is TRUE to
!     indicate that we are replacing time template strings only.
      fileName(2) = replaceTemplate( fileNameTemplate,.TRUE.,'find_file'  &
                                       ,fileDateName,fileTimeName )

!-------------------------------------------------------------------------------
    ELSE

!     NOT templateT

      IF ( nfile == 1 ) THEN
!       Assume that the required time lies in the first file.
        ifile = 1
      ELSE  !  nfile > 1

!       We have a list of file names and dates.

!       Initialise flag.
        found = .FALSE.
!-----------------------------------------------------------
        IF ( .NOT. climatol ) THEN

          DO ifile=1,nfile-1
            IF ( timeDate_cmp(time,date,'>='  &
                       ,fileTime(ifile),fileDate(ifile),'find_file')  .AND. &
                 timeDate_cmp(time,date,'<',  &
                       fileTime(ifile+1),fileDate(ifile+1),'find_file') ) THEN
              found = .TRUE.
              EXIT
           ENDIF
          ENDDO
!         If not yet found, check if time lies before start of first file.
          IF ( .NOT. found ) THEN
            IF ( timeDate_cmp( time,date,'<',fileTime(1),fileDate(1),'find_file' ) ) THEN
              ifile = 1   !   will calculate time relative to first file
            ELSE
!             We have not found a file covering this time, and know it is not before first file.
!             Assume that the required time lies in the last file (although we cannot check if
!             this is true without knowing end time of file).
              ifile = nfile
            ENDIF
          ENDIF   !  found

!-----------------------------------------------------------
        ELSE

!         climatol

!         Locate a file with the required date, ignoring the month.
!         If nfile=2 (meaning we have one file, plus a dummy file), then we always use
!         this first (climatological) file. Initialisation set date of dummy to be first date+1yr.
          IF ( nfile == 2 ) THEN
            ifile = 1
          ELSE
!           First remove year from date.
            CALL dateToBits( date,dateDay,dateMonth,dateYear,l_360,errMess )
            lookDate = dateMonth*100 + dateDay
            DO ifile=1,nfile-1
!             Remove year from file dates.
              CALL dateToBits( fileDate(ifile),dateDay,dateMonth,dateYear,l_360,errMess )
              tmpDate1 = dateMonth*100 + dateDay
              CALL dateToBits( fileDate(ifile+1),dateDay,dateMonth,dateYear,l_360,errMess )
              tmpDate2 = dateMonth*100 + dateDay
              IF ( lookDate>=tmpDate1 .AND. lookDate<tmpDate2 ) THEN
                found = .TRUE.
                EXIT
             ENDIF
            ENDDO
!           If haven't found a file, use last file.
            IF ( .NOT. found ) ifile = nfile
          ENDIF   !  nfile

!         Change dates of all files (other than last) to use required year.
!         Last file should use required year+1.
          CALL dateToBits( date,dateDay,dateMonth,nyr,l_360,errMess )
          DO i=1,nfile
            iyr = nyr
            IF ( i == nfile ) iyr = nyr + 1
            CALL dateToBits( fileDate(i),dateDay,dateMonth,dateYear,l_360,errMess )
            fileDate(i) = iyr*10000 + dateMonth*100 + dateDay
          ENDDO

        ENDIF  !  climatol
!-----------------------------------------------------------

      ENDIF !  nfile
!-------------------------------------------------------------------------------
    ENDIF  !  templateT
!-------------------------------------------------------------------------------
  ENDIF  !  next

!-------------------------------------------------------------------------------
! If necessary, calculate time level within file.
! This is always needed if NOT next, but may also be needed for next.

  IF ( findT ) THEN
   
!   Calculate time level within the file (or before first file).
    IF ( dataPer >= 0 ) THEN
!     Get time difference from start of file to given time.
      CALL timeDate_diff( fileTime(ifile),fileDate(ifile),time,date,l_360,'find_file'  &
                   ,nsec,nmin,nhr,ndy )
      nsec = ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec
!     Get time from start of file as a number of data intervals.
      t = FLOOR( REAL(nsec) / REAL(dataPer*itimeStep) ) + 1
   ELSEIF ( dataPer == periodMon ) THEN
!    Get time difference from start of file to given time.
     CALL timeDate_diff( fileTime(ifile),fileDate(ifile),time,date,l_360,'find_file'  &
                   ,nsec,nmin,nhr,ndy,nmon,nyr )
!    Convert to number of months (which equals number of data periods).
     t = FLOOR( REAL(nyr*12+nmon) ) + 1
   ENDIF
  ENDIF  !  findT

  END SUBROUTINE find_file

!###############################################################################
!###############################################################################
!###############################################################################
! subroutine varInfo
! Internal procedure in module misc_utils.
! Prints a message describing source of data for a variable.

  SUBROUTINE varInfo( varName,varFlag,varCode,varNameFile,fileFormat,constVal  &
                     ,readDump,varDesc,extraDesc )

  USE inout, ONLY :  &
!  imported scalar parameters
     formatNc,formatPP
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    varCode   &!  field code to use (e.g. STASH or PP code)
   ,varFlag    !  flag indicating source of data
!        >0 = field number
!         0 = variable not used (or a speical case)
!        -1 = use a constant value  <-1 = "special" cases

  INTEGER ::  &!  local SCALARS
    i,j   !  work

  REAL, INTENT(in), OPTIONAL ::  &!  optional in SCALARS
    constVal  !  a constant value

  LOGICAL, INTENT(in), OPTIONAL ::  &!  optional in SCALARS
    readDump  !  TRUE if field is to be read from a dump file (currently only
!                     used for initial condition)

  LOGICAL ::  &!  scalars
    fromDump  &!  work
   ,haveExtra  !  work

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    fileFormat   &!  format of file
   ,varName      &!  name of variable
   ,varNameFile   !  name of variable as used in a self-describing file

  CHARACTER(len=*), INTENT(in),OPTIONAL ::  &!  optional in SCALARS
    extraDesc  &!  more information (only used if varFlag<1)
   ,varDesc     !  a description of the variable

  CHARACTER(len=150) ::  &!  local SCALARS
    line      !  work
!-------------------------------------------------------------------------------

  fromDump = .FALSE.
  IF ( PRESENT(readDump) ) fromDump=readDump

  haveExtra = .FALSE.
  IF ( PRESENT(extraDesc) ) THEN
    IF ( LEN_TRIM(extraDesc) /= 0 ) haveExtra = .TRUE.
  ENDIF

! Echo variable's description.
  IF ( PRESENT(varDesc) ) WRITE(*,*) '# ',TRIM(varDesc)

! Load further information.
  line = ''

! Originally varFlag=0 was ignored, but now we echo further info if provided
! in extraDesc.

  IF ( varFlag/=0 .OR. (varFlag==0 .AND. haveExtra) ) line = TRIM(varName)

  IF ( varFlag > 0 ) THEN

    IF ( fromDump ) THEN
      line = TRIM(line) // ' is read from dump file.'
    ELSE

      IF ( fileFormat == formatPP ) THEN
        line =  TRIM(line) // ' uses field code='
        i = LEN_TRIM(line)
        WRITE(line(i:i+3),"(i4)") varCode
      ENDIF
      IF ( fileFormat /= formatNc ) THEN
        line = TRIM(line) // ' starts at field #'
        i = LEN_TRIM(line)
        WRITE(line(i:i+2),"(i3)") varFlag
      ELSE
        line = TRIM(line) // ' is called'
        i = LEN_TRIM(line)
        j = LEN_TRIM(varNameFile) + 1
        WRITE(line(i+1:i+j),"(a)") ' ' // TRIM(varNameFile)
        line = TRIM(line) // ' in input file'
      ENDIF

    ENDIF  !  fromDump

  ELSEIF ( varFlag == 0 ) THEN
    
    IF ( haveExtra ) line = TRIM(line) // ' ' // TRIM(extraDesc)

  ELSEIF ( varFlag == -1 ) THEN

    IF ( PRESENT(constVal) ) THEN
      line = TRIM(line) // ' is everywhere set to'
      i = LEN_TRIM(line)
      WRITE(line(i+1:i+9),"(es9.2)") constVal
    ELSE
      line = TRIM(line) // ' is everywhere set to a constant value'
    ENDIF

  ELSEIF ( varFlag < -1 ) THEN

    line = TRIM(line) // ' is set using flag='
    i = LEN_TRIM(line)
    WRITE(line(i+1:i+3),"(i2)") varFlag
    IF ( haveExtra ) line = TRIM(line) // ' ' // TRIM(extraDesc)

  ENDIF

  IF ( LEN_TRIM(line) > 0 ) WRITE(*,"(tr3,a)") TRIM(line)

  END SUBROUTINE varInfo
!###############################################################################
!###############################################################################
! subroutine varValue1D
! Summarise a 1-D field by printing range of values etc.

  SUBROUTINE varValue1D( summary,var,nval,varFormat,varDesc,varName )

  IMPLICIT NONE

! Scalar parameters.
  CHARACTER(len=12), PARAMETER :: defaultFm = '  tr1,es10.4'   !  default format for a single value

! Scalar arguments with intent(in)
  LOGICAL, INTENT(in) :: summary  !  flag indicating that summary statistics to
!         be printed.   T: print summary statistics (e.g. range)
!                       F: print values

! Array arguments with intent(in)
  REAL, INTENT(in) :: var(:)     !   the field to be described

! Optional scalar arguments with intent(in)
  INTEGER, INTENT(in), OPTIONAL :: nval   !  number of values to print per line
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varFormat !  format for write (for a
!                    single value of the field)
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varName   !  short name of field
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varDesc   !  description of field

! Local scalars.
  INTEGER :: nv   !  number of values to print per line
  INTEGER :: nx   !  size of variable
  CHARACTER(len=LEN(defaultFm)) :: valFm  !  format for a single value
  CHARACTER(len=3) :: cval   !  work
  CHARACTER(len=100) :: ff   !  a format

!-------------------------------------------------------------------------------
! Establish size of field.
!-------------------------------------------------------------------------------
  nx = SIZE( var,1 )

!-------------------------------------------------------------------------------
! Decide on format for a single value.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varFormat) ) THEN
    IF ( LEN_TRIM(varFormat) > LEN(valFm) ) THEN
!     The provided format is too long for available space, so use default format.
      valFm =defaultFm
    ELSE
      valFm = varFormat
    ENDIF
  ELSE
    valFm = defaultFm
  ENDIF

!-------------------------------------------------------------------------------
! Decide on number of values to print per line.
!-------------------------------------------------------------------------------
  IF ( PRESENT(nval) ) THEN
    nv = nval
  ELSE
    nv = 8
  ENDIF

!-------------------------------------------------------------------------------
! Print any available descriptors.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varName) ) THEN
    WRITE(*,"(2a)") '### name: ',TRIM(varName)
  ELSE
    WRITE(*,"(a)") '### name: not known'
  ENDIF
  IF ( PRESENT(varDesc) ) WRITE(*,"(a)") TRIM(varDesc)
  WRITE(*,"(a,i5)")'Size of field =',nx

!-------------------------------------------------------------------------------
! Print values.
!-------------------------------------------------------------------------------

  IF ( summary ) THEN

!-------------------------------------------------------------------------------
!   Print summary statistics.
!-------------------------------------------------------------------------------

    ff = "(a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
    WRITE(*,ff) 'Range of field=',MINVAL( var(:) ),' to ',MAXVAL( var(:) )

  ELSE

!-------------------------------------------------------------------------------
!   Print values of field.
!-------------------------------------------------------------------------------

!   Get a format such as ((10f6.1))
    WRITE(cval,"(i3)") nv
    ff = "((" // TRIM(cval) // TRIM(valFm) // "))"
    WRITE(*,ff) var(:)

  ENDIF  !  summary

  END SUBROUTINE varValue1D
  
!###############################################################################
!###############################################################################
! subroutine varValue2D
! Summarise a 2-D field by printing range of values etc.
! The 2nd dimension is expected to be levels (or pseudo-levels).

  SUBROUTINE varValue2D( doLevs,summary,var,nval,varFormat,varDesc,varName  &
                        ,maskMin,maskVar,levName )

  IMPLICIT NONE

! Scalar parameters.
  CHARACTER(len=12), PARAMETER :: defaultFm = '  tr1,es10.4'   !  default format for a single value

! Scalar arguments with intent(in)
  LOGICAL, INTENT(in) :: doLevs   !  flag indicating if each level is to be
!       treated separately   T: consider each level separately
!                            F: consider all levels together (i.e. whole field)
  LOGICAL, INTENT(in) :: summary  !  flag indicating that summary statistics to
!         be printed.   T: print summary statistics (e.g. range)
!                       F: print values

! Array arguments with intent(in)
  REAL, INTENT(in) :: var(:,:)   !   the field to be described

! Optional scalar arguments with intent(in)
  INTEGER, INTENT(in), OPTIONAL :: nval   !  number of values to print per line
  REAL, INTENT(in), OPTIONAL :: maskMin   !  minimum value when masking out
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varDesc   !  description of field
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varFormat !  format for write (for a
!                    single value of the field)
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varName   !  short name of field

! Optional array arguments with intent(in)
  REAL, INTENT(in), OPTIONAL :: maskVar(:,:)  !  field to use when masking out
!       Only values > maskMin are considered when summarising field.
!       Only used if summary=T AND doLevs=TRUE.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: levName(:)   !  name for each level
!                    e.g. name of PFT. Size should equal size(var,2).

! Local scalars.
  INTEGER :: iz   !  loop counter
  INTEGER :: n    !  work
  INTEGER :: nv   !  number of values to print per line
  INTEGER :: nx   !  size of variable in 1st dimension (expected to be "space"
!                       points)
  INTEGER :: nz   !  size of variable in 2nd dimension (expected to be "levels")
  REAL :: valMin  !  work
  REAL :: valMax  !  work
  LOGICAL :: doLevels  !  flag indicating if levels to be treated separately
  CHARACTER(len=LEN(defaultFm)) :: valFm  !  format for a single value
  CHARACTER(len=15) :: cval   !  work
  CHARACTER(len=30) :: clab   !  work
  CHARACTER(len=100) :: ff   !  a format

!-------------------------------------------------------------------------------
! Establish size of field.
!-------------------------------------------------------------------------------
  nx = SIZE( var,1 )
  nz = SIZE( var,2 )

  doLevels = doLevs
! For a single point field, don't treat levels separately.
  IF ( nx == 1 ) doLevels = .FALSE.

!-------------------------------------------------------------------------------
! Decide on format for a single value.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varFormat) ) THEN
    IF ( LEN_TRIM(varFormat) > LEN(valFm) ) THEN
!     The provided format is too long for available space, so use default format.
      valFm =defaultFm
    ELSE
      valFm = varFormat
    ENDIF
  ELSE
    valFm = defaultFm
  ENDIF

!-------------------------------------------------------------------------------
! Decide on number of values to print per line.
!-------------------------------------------------------------------------------
  IF ( PRESENT(nval) ) THEN
    nv = nval
  ELSE
    nv = 8
  ENDIF

!-------------------------------------------------------------------------------
! Print any available descriptors.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varName) ) THEN
    WRITE(*,"(2a)") '### name: ',TRIM(varName)
  ELSE
    WRITE(*,"(a)") '### name: not known'
  ENDIF
  IF ( PRESENT(varDesc) ) WRITE(*,"(a)") TRIM(varDesc)
  WRITE(*,"(2(a,i5))")'Size of field (x,z)=',nx,' x',nz

!-------------------------------------------------------------------------------
! Print values.
!-------------------------------------------------------------------------------

  IF ( summary ) THEN

!-------------------------------------------------------------------------------
!   Print summary statistics.
!-------------------------------------------------------------------------------

    IF ( doLevels ) THEN

!     Process each level separately.
      IF ( PRESENT(levName) ) THEN
        ff = "(a,i3,2a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
      ELSE
        ff = "(a,i3,a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
      ENDIF

      DO iz=1,nz

!       Get label.
        WRITE(clab,"(a7,i3)") 'Level #',iz
        IF ( PRESENT(levName) ) clab = TRIM(clab) // ' ' // TRIM(levName(iz))

!       Get range.
        IF ( PRESENT(maskVar) ) THEN
          n = COUNT( maskVar(:,iz) >= maskMin )
          IF ( n > 0 ) THEN
            valMin = MINVAL( var(:,iz), maskVar(:,iz)>=maskMin )
            valMax = MAXVAL( var(:,iz), maskVar(:,iz)>=maskMin )
          ENDIF
        ELSE
          n = 1
          valMin = MINVAL( var(:,iz) )
          valMax = MAXVAL( var(:,iz) )
        ENDIF

        IF ( n == 0 ) THEN
          WRITE(*,"(2a)") clab,': no values'
        ELSE
          WRITE(*,ff) clab,'Range =',valMin,' to ',valMax
        ENDIF

      ENDDO

    ELSE

!     Process all levels together.
      IF ( PRESENT(maskVar) ) THEN
        n = COUNT( maskVar(:,:) >= maskMin )
        IF ( n > 0 ) THEN
          valMin = MINVAL( var(:,:), maskVar(:,:)>=maskMin )
          valMax = MAXVAL( var(:,:), maskVar(:,:)>=maskMin )
        ENDIF
      ELSE
        n = 1
        valMin = MINVAL( var(:,:) )
        valMax = MAXVAL( var(:,:) )
      ENDIF

      IF ( n == 0 ) THEN
        WRITE(*,"(a)") 'no values'
      ELSE
        ff = "(a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
        WRITE(*,ff) 'Range of field =',valMin,' to ',valMax
      ENDIF

    ENDIF   !   levels

  ELSE

!-------------------------------------------------------------------------------
!   Print values of field.
!-------------------------------------------------------------------------------

!   Get a format such as ((10f6.1))
    WRITE(cval,"(i3)") nv
    ff = "((" // TRIM(cval) // TRIM(valFm) // "))"

    IF ( doLevels ) THEN
!     Print each level separately.
      DO iz=1,nz
        IF ( PRESENT(levName) ) THEN
          WRITE(*,"(a,i3,tr1,a)") 'Level #',iz,TRIM(levName(iz))
        ELSE
          WRITE(*,"(a,i3)") 'Level #',iz
        ENDIF
        WRITE(*,ff) var(:,iz)
      ENDDO
    ELSE
!     Print the full field as one.
      WRITE(*,ff) var(:,:)
    ENDIF

  ENDIF  !  summary

  END SUBROUTINE varValue2D
  
!###############################################################################
!###############################################################################
! subroutine varValue3D
! Summarise a 3-D field by printing range of values etc.
! The 2nd dimension is expected to be tiles,. The 3rd is expected to be levels
! and these are ignored - i.e. all levels on a tile are treated together.
! Expected use is for snow layer variables A(land_pts,ntiles,nsmax).

  SUBROUTINE varValue3D( doLevs,summary,var,nval,varFormat,varDesc,varName  &
                        ,maskMin,maskVar,levName )

  IMPLICIT NONE

! Scalar parameters.
  CHARACTER(len=12), PARAMETER :: defaultFm = '  tr1,es10.4'   !  default format for a single value

! Scalar arguments with intent(in)
  LOGICAL, INTENT(in) :: doLevs   !  flag indicating if each level is to be
!       treated separately   T: consider each level separately
!                            F: consider all levels together (i.e. whole field)
!       In this context, levels are tiles.
  LOGICAL, INTENT(in) :: summary  !  flag indicating that summary statistics to
!         be printed.   T: print summary statistics (e.g. range)
!                       F: print values

! Array arguments with intent(in)
  REAL, INTENT(in) :: var(:,:,:)   !   the field to be described

! Optional scalar arguments with intent(in)
  INTEGER, INTENT(in), OPTIONAL :: nval   !  number of values to print per line
  REAL, INTENT(in), OPTIONAL :: maskMin   !  minimum value when masking out
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varDesc   !  description of field
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varFormat !  format for write (for a
!                    single value of the field)
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varName   !  short name of field

! Optional array arguments with intent(in)
  REAL, INTENT(in), OPTIONAL :: maskVar(:,:)  !  field to use when masking out
!       Only values > maskMin are considered when summarising field.
!       Only used if summary=T AND doLevs=TRUE.
  CHARACTER(len=*), INTENT(in), OPTIONAL :: levName(:)   !  name for each level
!                    e.g. name of PFT. Size should equal size(var,2).

! Local scalars.
  INTEGER :: iz   !  loop counter
  INTEGER :: j    !  work
  INTEGER :: n    !  work
  INTEGER :: nn   !  size of variable in 3rd dimension (expected to be snow layers)
  INTEGER :: nv   !  number of values to print per line
  INTEGER :: nx   !  size of variable in 1st dimension (expected to be "space"
!                       points)
  INTEGER :: nz   !  size of variable in 2nd dimension (expected to be "levels")
  REAL :: val     !  work
  REAL :: valMin  !  work
  REAL :: valMax  !  work
  LOGICAL :: doLevels  !  flag indicating if levels to be treated separately
  CHARACTER(len=LEN(defaultFm)) :: valFm  !  format for a single value
  CHARACTER(len=15) :: cval   !  work
  CHARACTER(len=30) :: clab   !  work
  CHARACTER(len=100) :: ff   !  a format

! Local arrays.
  REAL :: varSlice(SIZE(var,1),1)  !  work

!-------------------------------------------------------------------------------
! Establish size of field.
!-------------------------------------------------------------------------------
  nx = SIZE( var,1 )
  nz = SIZE( var,2 )
  nn = SIZE( var,3 )

  doLevels = doLevs
! For a single point and single level field, don't treat levels (dims #2 and 3)
! separately.
  IF ( nx*nz == 1 ) doLevels = .FALSE.

!-------------------------------------------------------------------------------
! Decide on format for a single value.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varFormat) ) THEN
    IF ( LEN_TRIM(varFormat) > LEN(valFm) ) THEN
!     The provided format is too long for available space, so use default format.
      valFm =defaultFm
    ELSE
      valFm = varFormat
    ENDIF
  ELSE
    valFm = defaultFm
  ENDIF

!-------------------------------------------------------------------------------
! Decide on number of values to print per line.
!-------------------------------------------------------------------------------
  IF ( PRESENT(nval) ) THEN
    nv = nval
  ELSE
    nv = 9
  ENDIF

!-------------------------------------------------------------------------------
! Print any available descriptors.
!-------------------------------------------------------------------------------
  IF ( PRESENT(varName) ) THEN
    WRITE(*,"(2a)") '### name: ',TRIM(varName)
  ELSE
    WRITE(*,"(a)") '### name: not known'
  ENDIF
  IF ( PRESENT(varDesc) ) WRITE(*,"(a)") TRIM(varDesc)
  WRITE(*,"(3(a,i5))")'Size of field=',nx,' x',nz,' x',nn

!-------------------------------------------------------------------------------
! Print values.
!-------------------------------------------------------------------------------

  IF ( summary ) THEN

!-------------------------------------------------------------------------------
!   Print summary statistics.
!-------------------------------------------------------------------------------

    IF ( doLevels ) THEN

!     Process each level separately.
      IF ( PRESENT(levName) ) THEN
        ff = "(a,i3,2a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
      ELSE
        ff = "(a,i3,a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
      ENDIF

      DO iz=1,nz

!       Get label.
        WRITE(clab,"(a7,i3)") 'Level #',iz
        WRITE(*,*)'1st clab=',TRIM(clab)
        IF ( PRESENT(levName) ) clab = TRIM(clab) // ' ' // TRIM(levName(iz))
        WRITE(*,*)'2nd clab=',TRIM(clab)

!       Get range.
        IF ( PRESENT(maskVar) ) THEN
          n = COUNT( maskVar(:,iz) >= maskMin )
          IF ( n > 0 ) THEN
            valMin = 1.0e20
            valMax = -1.0e20
!           Loop over 3rd dimension (e.g. snow layers) for this tile.
            DO j=1,nn
!              varSlice(:,1) =  reshape( var(:,iz,j), (/ nx /) )
              varSlice(:,1) =  var(:,iz,j)
              val = MINVAL( varSlice(:,1), maskVar(:,iz)>=maskMin )
              IF ( val < valMin ) valMin = val
              val = MAXVAL( varSlice(:,1), maskVar(:,iz)>=maskMin )
              IF ( val > valMax ) valMax = val
            ENDDO
          ENDIF
        ELSE
          n = 1
          valMin = MINVAL( var(:,iz,:) )
          valMax = MAXVAL( var(:,iz,:) )
        ENDIF

        IF ( n == 0 ) THEN
          WRITE(*,"(2a)") clab,': no values'
        ELSE
          WRITE(*,ff) clab,'Range =',valMin,' to ',valMax
        ENDIF

      ENDDO

    ELSE

!     Process all levels together.
      IF ( PRESENT(maskVar) ) THEN
        n = COUNT( maskVar(:,:) >= maskMin )
        IF ( n > 0 ) THEN
          valMin = 1.0e20
          valMax = -1.0e20
          DO iz=1,nz
            DO j=1,nn
              varSlice(:,1) =  var(:,iz,j)
              val = MINVAL( varSlice(:,1), maskVar(:,iz)>=maskMin )
              IF ( val < valMin ) valMin = val
              val = MAXVAL( varSlice(:,1), maskVar(:,iz)>=maskMin )
              IF ( val > valMax ) valMax = val
            ENDDO
          ENDDO
        ENDIF
      ELSE
        n = 1
        valMin = MINVAL( var(:,:,:) )
        valMax = MAXVAL( var(:,:,:) )
      ENDIF

      IF ( n == 0 ) THEN
        WRITE(*,"(a)") 'no values'
      ELSE
        ff = "(a," // TRIM(valFm) // ",a," // TRIM(valFm) // ")"
        WRITE(*,ff) 'Range of field =',valMin,' to ',valMax
      ENDIF

    ENDIF   !   levels

  ELSE

!-------------------------------------------------------------------------------
!   Print values of field.
!-------------------------------------------------------------------------------

!   Get a format such as ((10f6.1))
    WRITE(cval,"(i3)") nv
    ff = "((" // TRIM(cval) // TRIM(valFm) // "))"

    IF ( doLevels ) THEN
!     Print each level separately.
      DO iz=1,nz
        IF ( PRESENT(levName) ) THEN
          WRITE(*,"(a,i3,tr1,a)") 'Level #',iz,TRIM(levName(iz))
        ELSE
          WRITE(*,"(a,i3)") 'Level #',iz
        ENDIF
        DO n=1,nn
          WRITE(*,ff) var(:,iz,n)
        ENDDO
      ENDDO
    ELSE
!     Print the full field as one. Note that this won't appear in the most
!     intuitive order....but it's good enough.
      WRITE(*,ff) var(:,:,:)
    ENDIF

  ENDIF  !  summary

  END SUBROUTINE varValue3D
  
!###############################################################################
!###############################################################################
! subroutine check_template
! Internal procedure in module misc_utils.
! Given a file name, check if it appears to be a valid template name.
! Note that this has not been properly checked since MetO filenames accommodated
! elsewhere. In particular, decades are scarcely considered.

  SUBROUTINE check_template( dataPer,filePer,itimeStep,time,date,climatol  &
                            ,fullCheck,fileName,errMess  &
                            ,templateT,templateV,templateUnits )

  USE inout, ONLY :  &
!  imported scalar parameters
     ntemplString,periodAnn,periodMon,periodOneFile  &
!  imported array parameters
    ,templString

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_360

  USE time_mod, ONLY :  &
!  imported procedures
     dateToBits,get_period,s_to_chhmmss
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalars with intent(in)

  INTEGER, INTENT(in) ::  &!  in scalars
    filePer  &!  period of file
!                  >0 is period as number of timesteps
!                 <=0 are special cases
   ,dataPer  &!  period of data in file
!                  >0 is period in seconds
!                 <=0 are special cases
   ,itimeStep &!  timestep length (s) for filePer or dataPer>0
   ,time     &!  time of day (s) of first data in file
   ,date      !  date (yyyymmdd) of first data in file

  LOGICAL, INTENT(in) ::  &!  in scalars
    climatol   &!  T means files are to be treated as climatological (year unimportant)
   ,fullCheck   !  T means all aspects of the template are checked
!               !  F means a limited check (just to establish if templating is indicated)

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    fileName  &!  name of file (possibly including directory and template parts)
   ,errMess    !  message printed on error

!-------------------------------------------------------------------------------
! Scalars with intent(out)

  LOGICAL, INTENT(out) ::  &!  out scalars
    templateT   &!  T if file name is identified as suitable for time templating
   ,templateV    !  T if file name is identified as suitable for variable
!                      name templating

  CHARACTER(len=2), INTENT(out) ::  &!  out scalars
    templateUnits   !  indicates shortest time period that is included in a template

!-------------------------------------------------------------------------------
! Local scalars

  INTEGER ::  &!  local scalars
    dateDay,dateMonth,dateYear  &!  parts of date
   ,i,i1,i2,ILEN,ipos,j   !  work

  LOGICAL ::  &!  local scalars
    err1,err2,found     &!  work
   ,haveDay,haveHour,haveMin,haveMonth,haveYear  ! TRUE if template uses the named time unit

  CHARACTER(len=2) ::  &!  local scalars
    dataPerUnits  &!  units used to describe file period
   ,filePerUnits   !  units used to describe file period

  CHARACTER(len=LEN(templString)) ::  &!  local scalars
    tString    !  used to hold a template string

  CHARACTER(len=10) ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
!                                 ! that is required because some compilers can not cope
!                                 ! with recursive write statements
!-------------------------------------------------------------------------------

! Start by assuming that there are no templating strings.
  templateT = .FALSE.
  templateV = .FALSE.

! Look for '%' in the name.
  ipos = INDEX ( fileName, '%' )

! If there is no '%' in the string, there's no templating and no more to do.
  IF ( ipos == 0 ) RETURN

!-------------------------------------------------------------------------------
! Templating appears to be used.

! Initialise.
  templateUnits = ''
  haveYear = .FALSE.
  haveMonth = .FALSE.
  haveDay = .FALSE.
  haveHour = .FALSE.
  haveMin = .FALSE.

! Go through the string, searching for template strings.
! Note that these are all assumed to be the same length.
  i1 = ipos
  ILEN = LEN_TRIM( fileName )
  DO
    IF ( i1 > ILEN ) EXIT

!   Search for a %.
    ipos = INDEX( fileName(i1:ILEN), '%' )
    IF ( ipos == 0 ) EXIT
    ipos = ipos + i1 - 1   !  location in fileName
    i2 = MIN( ipos+LEN(templString)-1,ILEN )
    tString = fileName(ipos:i2)

!   Check we have a recognised template string.
    found = .FALSE.
    DO i=1,ntemplString

      IF ( tString == templString(i) ) THEN
!       Template string recognised.
        found = .TRUE.

!       Decide if this is time templating or not.
        IF ( templString(i) == '%vv' ) THEN
          templateV = .TRUE.
        ELSE
          templateT = .TRUE.
!         Work out what time unit this represents.
          SELECT CASE ( tString )
            CASE ( '%tc' )
!             Decades. Nothing more to do here.
              haveYear = .TRUE.
            CASE ( '%y2', '%y4', '%yc' )
              haveYear = .TRUE.
            CASE ( '%m1', '%m2', '%mc', '%mm' )
              haveMonth = .TRUE.
            CASE ( '%d1', '%d2', '%dc' )
              haveDay = .TRUE.
            CASE ( '%h1', '%h2', '%hc' )
              haveHour = .TRUE.
            CASE ( '%n2' )
              haveMin = .TRUE.
            CASE default
              WRITE(*,*)'ERROR: check_template: no code for template=',tString
              STOP
          END SELECT
        ENDIF  !  templString
        EXIT  !  loop over templString
      ENDIF  !  tstring

    ENDDO  !  templString

    IF ( .NOT. found ) THEN
      WRITE(*,*)'ERROR: check_template: do not recognise template string: ',tString
      WRITE(*,*)'fileName=',TRIM(fileName)
      WRITE(*,"(a,(15(tr1,a)))")'Valid template strings:',(TRIM(templString(j)),j=1,ntemplString)
      WRITE(*,*)TRIM( errMess )
      STOP
    ENDIF

!   Set index for next loop.
    i1 = ipos + LEN(templString)

  ENDDO   !  loop through string
!-------------------------------------------------------------------------------

! If a basic check was requested, we are finished.
  IF ( .NOT. fullCheck ) RETURN

!-------------------------------------------------------------------------------

! More detailed checking of template.

  IF ( templateT ) THEN

!   Climatological files must not include the year in time-templating.
    IF ( climatol .AND. haveYear ) THEN
      WRITE(*,*)'ERROR: check_template: climatological files must not include year'
      WRITE(*,*)'in time-templating.'
      WRITE(*,*)'fileName=',TRIM(fileName)
      WRITE(*,*)TRIM( errMess )
      STOP
    ENDIF
  
!   Get shortest period that appears in any template string.
    IF ( haveMin ) THEN
      templateUnits = 'mn'
    ELSEIF ( haveHour ) THEN
      templateUnits = 'hr'
    ELSEIF ( haveDay ) THEN
      templateUnits = 'dy'
    ELSEIF ( haveMonth ) THEN
      templateUnits = 'mo'
    ELSE
      templateUnits = 'yr'
    ENDIF
  
!   Work out what units to use for file period.
    CALL get_period( dataPer,itimeStep,i,dataPerUnits )
    IF ( filePer /= periodOneFile ) THEN
      CALL get_period( filePer,itimeStep,i,filePerUnits )
    ELSE
!     This period is unsuitable, and error is raised below.
      filePerUnits = 'XX'
    ENDIF
  
!   Check that units are acceptable for use with a template.
!   At time of writing, this was guaranteed, but check anyway.
    SELECT CASE ( fileperUnits )
      CASE ( 'yr', 'mo', 'dy', 'hr', 'mn' )
      CASE default
        WRITE(*,*)'ERROR: check_template: period of file not suitable for template.'
        WRITE(*,*)'fileName=',TRIM(fileName)
        IF ( filePer /= periodOneFile ) THEN
          WRITE(*,*)'Time units not available via template.'
          WRITE(*,*)'filePer=',filePer,' units=fileperUnits=',filePerUnits
        ELSE
          WRITE(*,*)'Cannot use time-templating with file period=periodOneFile=',periodOneFile
          WRITE(*,*)'Change file period, or don''t use time templating.'
        ENDIF
        WRITE(*,*) TRIM( errMess )
        STOP
    END SELECT
  
!   Check that the units identified for file period appear in the template.
!   Also check that shorter time periods are not found in the template.
    err1 = .FALSE.
    err2 = .FALSE.
    SELECT CASE ( filePerUnits )
      CASE ( 'yr' )
        IF ( .NOT. haveYear ) err1 = .TRUE.
        IF ( haveMonth .OR. haveDay .OR. haveHour .OR. haveMin ) err2 = .TRUE.
      CASE ( 'mo' )
        IF ( .NOT. haveMonth ) err1 = .TRUE.
        IF ( haveDay .OR. haveHour .OR. haveMin ) err2 = .TRUE.
      CASE ( 'dy' )
        IF ( .NOT. haveDay ) err1 = .TRUE.
        IF ( haveHour .OR. haveMin ) err2 = .TRUE.
      CASE ( 'hr' )
        IF ( .NOT. haveHour ) err1 = .TRUE.
        IF ( haveMin ) err2 = .TRUE.
      CASE ( 'mn' )
        IF ( .NOT. haveMin ) err1 = .TRUE.
      CASE default
        WRITE(*,*)'ERROR: check_template: no code for filePerUnits=',filePerUnits
        WRITE(*,*)'fileName=',TRIM(fileName)
        WRITE(*,*) TRIM( errMess )
        STOP
    END SELECT
  
    IF ( err1 ) THEN
      WRITE(*,*)'ERROR: check_template: units for file period not in template.'
      WRITE(*,*)'fileName=',TRIM(fileName)
      WRITE(*,*)'fileperUnits=',filePerUnits
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF
    IF ( err2 ) THEN
      WRITE(*,*)'ERROR: check_template: shorter time periods in template.'
      WRITE(*,*)'filePerUnits=',filePerUnits,' but we have found shorter time units.'
      WRITE(*,*)'fileName=',TRIM(fileName)
      WRITE(*,*)'haveYear,month,day,hour,min=',haveYear,haveMonth,haveDay,haveHour,haveMin
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF
  
!-------------------------------------------------------------------------------
!   Check that the time and date associated with the template is acceptable.
!   This is awkward - haven't worked out a simpler solution yet.
!   This is coded assuming that dataPer is either >0, or is one of periodMon
!   or periodAnn, so test that.
    IF ( dataPer<=0 .AND. dataPer/=periodMon .AND. dataPer/=periodAnn ) THEN
      WRITE(*,*)'ERROR: check_template: dataPer not expected! See code.'
      STOP
    ENDIF
!   Split date into components.
    CALL dateToBits( date,dateDay,dateMonth,dateYear,l_360,errMess )
    err1 = .FALSE.
    SELECT CASE ( filePerUnits )
  
!-----------------------------------------------------------------
      CASE ( 'yr' )
        SELECT CASE ( dataPerUnits )
          CASE ( 'yr' )
!           No restriction.
          CASE ( 'mo' )
!           Date should be in January.
            IF ( dateMonth /= 1 ) THEN
              WRITE(*,*)'Date should be in January.'
              err1 = .TRUE.
            ENDIF
!           For simplicity, want a day that is in all months - e.g. not 31.)
           IF ( .NOT. l_360 ) THEN
!            If day of month > 28, this can be difficult - so avoid.
             IF ( dateDay > 28 ) THEN
               WRITE(*,*)'ERROR: check_template: dateDay > 28.'
               WRITE(*,*)'The date given has day of month > 28.'
               WRITE(*,*)'This is difficult, because not all months have this number of days!'
               WRITE(*,*)'More code is probably required!'
               WRITE(*,*)'If data are timestamped for end of month, use 00H on 1st of next month.'
               STOP
             ENDIF
           ELSE   !  l_360
             IF ( dateDay > 30 ) THEN
              WRITE(*,*)'ERROR: check_template: dateDay > 30.'
              WRITE(*,*)'The date given has day of month > 30.'
              WRITE(*,*)'This cannot be used with 360-day calendar.'
              STOP
            ENDIF
          ENDIF  !  l_360
          CASE ( 'dy' )
!           Date should be 01Jan.
            IF ( dateMonth/=1 .OR. dateDay/=1 ) THEN
              WRITE(*,*)'Date should be 1st January.'
              err1 = .TRUE.
            ENDIF
          CASE ( 'hr', 'mn' )
!           Date should be 00H 01Jan.
            IF ( dateMonth/=1 .OR. dateDay/=1 .OR. time/=0 ) THEN
              WRITE(*,*)'Date should be 00H 1st January.'
              err1 = .TRUE.
            ENDIF
        END SELECT
!  ---------------------------------------------------
      CASE ( 'mo' )
        SELECT CASE ( dataPerUnits )
          CASE ( 'mo' )
!           No restriction.
          CASE ( 'dy' )
!           Date should be 1st of a month.
            IF ( dateDay/=1 ) THEN
              WRITE(*,*)'Date should be 1st of a month.'
              err1 = .TRUE.
            ENDIF
          CASE ( 'hr', 'mn' )
!           Date should be 00H on 1st of a month.
            IF ( dateDay/=1 .OR. time/=0 ) THEN
              WRITE(*,*)'Date should be 00H on 1st of a month.'
              err1 = .TRUE.
            ENDIF
        END SELECT
!  ---------------------------------------------------
      CASE ( 'dy' )
        SELECT CASE ( dataPerUnits )
          CASE ( 'dy' )
!           No restriction.
          CASE ( 'hr', 'mn' )
!           Time should be 00H.
            IF ( time /= 0 ) THEN
              WRITE(*,*)'Time should be 00H.'
              err1 = .TRUE.
            ENDIF
        END SELECT
!  ---------------------------------------------------
      CASE ( 'hr' )
        SELECT CASE ( dataPerUnits )
          CASE ( 'hr' )
!           No restriction.
          CASE ( 'mn' )
!           Time should be 00H.
            IF ( time /= 0 ) THEN
              WRITE(*,*)'Time should be 00H.'
              err1 = .TRUE.
            ENDIF
        END SELECT
!  ---------------------------------------------------
      CASE ( 'mn' )
!       dataPerUnits must be 'mn', and there is no restriction on time or date.
!  ---------------------------------------------------
      CASE default
        WRITE(*,*)'ERROR: check_template: no code for filePerUnits=',filePerUnits
        WRITE(*,*)'fileName=',TRIM(fileName)
        WRITE(*,*) TRIM( errMess )
        STOP
    END SELECT
!  -------------------------------------------------------------------------------
  
    IF ( err1 ) THEN
      WRITE(*,*)'ERROR: check_template: time or date of file is not suitable.'
      WRITE(*,*)'For template files, the time and date given for the files needs'
      WRITE(*,*)'to be consistent with the periods of the data and the files.'
      WRITE(*,*)'fileName=',TRIM(fileName)
      time_hms = s_to_chhmmss( time )
      WRITE(*,*)'time=',time,'s (',time_hms,'), date=',date
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF

  ENDIF  !  templateT

  END SUBROUTINE check_template

!###############################################################################
!###############################################################################
!#######################################################################
!#######################################################################
!#######################################################################
! subroutine read_list
! Read lines from standard input and parse into lists.

  SUBROUTINE read_list( inUnit,nfMax,nlMax,startTag,endTag,delim,errMess,nline  &
                       ,ivar,ivarPos,rvar,rvarPos,cvar1,cvar1Pos,cvar2,cvar2Pos &
                       ,cvar3,cvar3Pos,cvar4,cvar4Pos  &
                       ,lvar1,lvar1Pos,lvar2,lvar2Pos )

! Read lines and parse into lists.
! Lines are known to lie between two known "tag" lines.
! Tailored for current needs of JULES!
! Use keywords for optional arguments.
!-------------------------------------------------------------------------------

  USE file_utils, ONLY :  &
!  imported procedures
     findTag

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar parameters.
  INTEGER, PARAMETER ::  &
    clen = 150   !  length of character variable (needs to be large enough for
!                                                 longest field expected)

!-------------------------------------------------------------------------------
! imported scalars with intent(in)

  INTEGER, INTENT(in) ::  &
    inUnit  &!  unit to read from
   ,nfMax   &!  maximum number of fields to read per line
   ,nlMax    !  maximum possible number of lines to read

  INTEGER, OPTIONAL, INTENT(in) :: &
    ivarPos      &!  field number associated with ivar
   ,cvar1Pos     &!  field number associated with cvar1
   ,cvar2Pos     &!  field number associated with cvar2
   ,cvar3Pos     &!  field number associated with cvar3
   ,cvar4Pos     &!  field number associated with cvar4
   ,lvar1Pos     &!  field number associated with lvar1
   ,lvar2Pos     &!  field number associated with lvar2
   ,rvarPos       !  field number associated with rvar

  CHARACTER(len=1), INTENT(in) ::  &
    delim     !  character used to delimit (separate) values in list

  CHARACTER(len=*), INTENT(in) ::  &
    startTag &!  tag indicating start of list
   ,endTag   &!  tag indicating end of list
   ,errMess   !  message printed if error raised

!-------------------------------------------------------------------------------
! imported scalars with intent(out)
  INTEGER, INTENT(out) ::  &
    nline   !  number of lines (not empty) read

!-------------------------------------------------------------------------------
! imported optional arrays with intent(out)

  INTEGER, OPTIONAL, INTENT(out) ::  &
    ivar(nlMax)  !  an integer variable to be found in line

  REAL, OPTIONAL, INTENT(out) ::  &
    rvar(nlMax)  !  a real variable to be found in line

  LOGICAL, OPTIONAL, INTENT(out) :: lvar1(nlMax)  !  a logical variable to be found in line
  LOGICAL, OPTIONAL, INTENT(out) :: lvar2(nlMax)  !  a logical variable to be found in line

  CHARACTER(len=*), OPTIONAL, INTENT(out)  ::  &
    cvar1(nlMax)  &!  a character variable to be found in line
   ,cvar2(nlMax)  &!  a character variable to be found in line
   ,cvar3(nlMax)  &!  a character variable to be found in line
   ,cvar4(nlMax)   !  a character variable to be found in line
!-------------------------------------------------------------------------------
! local scalars

  INTEGER ::  &
    ierr   &!  erorr value
   ,wlen    !  work

  LOGICAL ::  &
    errFlag  !  TRUE indicates an error has been raised

  CHARACTER(len=20) ::  &
    cFormat  !  used to hold a format

  CHARACTER(len=200) ::  &
    inLine   !  a line read from file.  
!------------------------------------------------------------------------------- 
! local scalars

  INTEGER ::  &
    nfieldFound &!  number of fields found by getFields
   ,nvar         !  number of variables requested
!-------------------------------------------------------------------------------
! local arrays

  CHARACTER(len=clen) ::  &
    fieldList(nfMax)   !  fields in line

!-------------------------------------------------------------------------------

! Work out how many variables have been passed (and are to be filled).
  nvar = 0
  IF ( PRESENT( iVar ) ) nvar = nvar + 1
  IF ( PRESENT( rVar ) ) nvar = nvar + 1
  IF ( PRESENT( cVar1 ) ) nvar = nvar + 1
  IF ( PRESENT( cVar2 ) ) nvar = nvar + 1
  IF ( PRESENT( cVar3 ) ) nvar = nvar + 1
  IF ( PRESENT( cVar4 ) ) nvar = nvar + 1
  IF ( PRESENT( lVar1 ) ) nvar = nvar + 1
  IF ( PRESENT( lVar2 ) ) nvar = nvar + 1

! Check that optional arguments are provided in pairs.
  IF ( PRESENT(ivar) .NEQV. PRESENT(ivarPos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'ivar and ivarPos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(rvar) .NEQV. PRESENT(rvarPos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'rvar and rvarPos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(cvar1) .NEQV. PRESENT(cvar1Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'cva1r and cvar1Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(cvar2) .NEQV. PRESENT(cvar2Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'cvar2 and cvar2Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(cvar3) .NEQV. PRESENT(cvar3Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'cvar3 and cvar3Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(cvar4) .NEQV. PRESENT(cvar4Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'cvar4 and cvar4Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(lvar1) .NEQV. PRESENT(lvar1Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'lvar1 and lvar1Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  IF ( PRESENT(lvar2) .NEQV. PRESENT(lvar2Pos) ) THEN
    WRITE(*,*)'ERROR: read_list: optional arguments not in pairs.'
    WRITE(*,*)'lvar2 and lvar2Pos should BOTH be provided or NEITHER should be provided.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
  
! Check that the numbers of the requested fields lie within range.
! We do not check for repeated values.
  errFlag = .FALSE.
  IF ( PRESENT(ivarPos)) THEN
     IF (ivarPos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: ivarPos>nfMax: ivarPos=',ivarPos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(rvarPos)) THEN
     IF (rvarPos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: rvarPos>nfMax: rvarPos=',rvarPos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(cvar1Pos)) THEN
     IF (cvar1Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: cvar1Pos>nfMax: cvar1Pos=',cvar1Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(cvar2Pos)) THEN
     IF (cvar2Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: cvar2Pos>nfMax: cvar2Pos=',cvar2Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(cvar3Pos)) THEN
     IF (cvar3Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: cvar3Pos>nfMax: cvar3Pos=',cvar3Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(cvar4Pos)) THEN
     IF (cvar4Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: cvar4Pos>nfMax: cvar4Pos=',cvar4Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(lvar1Pos)) THEN
     IF (lvar1Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: lvar1Pos>nfMax: lvar1Pos=',lvar1Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( PRESENT(lvar2Pos)) THEN
     IF (lvar2Pos>nfMax ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: lvar2Pos>nfMax: lvar2Pos=',lvar2Pos,' nfMax=',nfMax
     ENDIF
  ENDIF
  IF ( errFlag ) THEN
    WRITE(*,*)'ERROR: read_list: field location out of range.'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF
!-------------------------------------------------------------------------------
  nline = 0

! Find the tag that indicates the start of the list.
  CALL findTag( inUnit,errMess,startTag,preInit=.TRUE. )

!-------------------------------------------------------------------------------
! Read lines, stopping when end tag is found.
  DO

    READ(inUnit,fmt="(a)",iostat=ierr) inLine
    IF ( ierr /= 0 ) THEN
      WRITE(*,*)'ERROR: read_list: error on read.'
      IF ( ierr < 0 ) WRITE(*,*)'Error: end of file (iostat=',ierr,')'
      WRITE(*,*) TRIM( errMess ) 
      STOP
    ENDIF

!   Do nothing with empty lines.
    IF ( LEN_TRIM(inLine) == 0 ) CYCLE

!   Parse into nfMax or fewer fields.
    CALL getFields( clen,.TRUE.,.TRUE.,.FALSE.,delim,inLine,errFlag,nfMaxUse=.TRUE.  &
                   ,nfield=nfieldFound,fieldList=fieldList )

    IF ( errFlag ) THEN
      WRITE(*,*)'ERROR: read_list: error raised in getFields.'
      WRITE(*,*)'inLine=',TRIM(inLine)
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF

!   Exit when end tag is found.
    IF ( fieldList(1) == endTag ) EXIT

!   Check all required fields were found.
    IF ( nfieldFound < nvar ) THEN
      WRITE(*,*)'ERROR: read_list: not enough fields'
      WRITE(*,*)'Looking for ',nvar,' fields, found ',nfieldFound
      WRITE(*,*)'inLine=',TRIM(inLine)
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF

!   Increment counter.
    nline = nline + 1

    IF ( nline > nlMax ) THEN
      WRITE(*,*)'ERROR: read_list: nline > nlineMax.'
      WRITE(*,*)'List is too long for available storage.'
      WRITE(*,*)'Have exceeded max possible number of lines that can be read=',nlMax
      WRITE(*,*)'without finding end tag=',TRIM(endTag)
      WRITE(*,*) TRIM( errMess )
      STOP
    ENDIF

!   Parse into variables.
    IF ( PRESENT( ivar ) ) THEN
      wlen = LEN_TRIM( fieldList(ivarPos) )
      WRITE(cformat,"('(i',i1,')')") wlen
      READ(fieldList(ivarPos),cformat) ivar(nline)
    ENDIF
    IF ( PRESENT( rvar ) ) THEN
      wlen = LEN_TRIM( fieldList(rvarPos) )
      WRITE(cformat,"('(f',i2,'.0)')") wlen
      READ(fieldList(rvarPos),cformat) rvar(nline)
    ENDIF
    IF ( PRESENT( cvar1 ) ) cvar1(nline) = fieldList(cvar1Pos)
    IF ( PRESENT( cvar2 ) ) cvar2(nline) = fieldList(cvar2Pos)
    IF ( PRESENT( cvar3 ) ) cvar3(nline) = fieldList(cvar3Pos)
    IF ( PRESENT( cvar4 ) ) cvar4(nline) = fieldList(cvar4Pos)
    IF ( PRESENT( lvar1 ) ) READ(fieldList(lvar1Pos),"(l1)") lvar1(nline)
    IF ( PRESENT( lvar2 ) ) READ(fieldList(lvar2Pos),"(l1)") lvar2(nline)
      
  ENDDO


  END SUBROUTINE read_list
!###############################################################################
!###############################################################################
! subroutine get_filePer_vars
! Internal procedure in module misc_utils.
! Given file period flag, returns units and increment.

  SUBROUTINE get_filePer_vars( filePer,itimeStep,inc,units )

  USE inout, ONLY :  &!
!   imported scalar parameters
     periodAnn,periodMon

  IMPLICIT NONE
!-------------------------------------------------------------------------------

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) ::  &
    filePer    &!  file period
   ,itimeStep   !  model timestep (s)
!                  >0 indicates period as a number of timesteps
!                 <=0 indicate "special" cases (e.g. monthly).

! Scalar arguments with intent(out)
  CHARACTER(len=*), INTENT(out) ::  &
    units   !  units of time used for filePer

  INTEGER, INTENT(out) ::  &
    inc   !   increment represented by filePer
!-------------------------------------------------------------------------------
  IF ( filePer >= 0 ) THEN
    units = 'sec'
    inc = filePer*itimeStep
  ELSE
!   Get units and increment for special cases.
!   Currently these all use 'months' for units of increment.
    units = 'mon'
    SELECT CASE ( filePer )
      CASE ( periodAnn )
        inc = 12
      CASE ( periodMon )
        inc = 1
      CASE default
        WRITE(*,*)'ERROR: get_filePer_vars: no code for filePer=',filePer
      STOP
    END SELECT
  ENDIF  !  filePer

  END SUBROUTINE get_filePer_vars

!###############################################################################
!###############################################################################

  SUBROUTINE varList( nvarIn,sdfFile,varNameIn,varNameList  &
                     ,errFound,foundVar  &
                     ,varFlagIn,varNzIn,varConstIn,varLog1In,varLog2In  &
                     ,varFileNameIn,varInterpIn,varNameSDFin  &
                     ,varFlag,varNz,varConst,varLog1,varLog2  &
                     ,varFileName,varInterp,varNameSDF )

! A list of variable names (expected to have been read from an input file) is
! compared with a list of possible names in a "master" list. When a match is
! found, other values are saved in the order of the master list.
! All very hardwired to current needs.

!------------------------------------------------------------------------------- 
! Scalar arguments with intent(in)
!------------------------------------------------------------------------------- 
  INTEGER, INTENT(in) :: nvarIn  !  number of variables to process
  LOGICAL, INTENT(in) :: sdfFile !  TRUE if dealing with an SDF (self-describing file)

!------------------------------------------------------------------------------- 
! Array arguments with intent(in)
!------------------------------------------------------------------------------- 
  CHARACTER(len=*), INTENT(in) :: varNameIn(:)    !  list of variables names
  CHARACTER(len=*), INTENT(in) :: varNameList(:)  !  master list of (all possible) names

!------------------------------------------------------------------------------- 
! Scalar arguments with intent(out)
!------------------------------------------------------------------------------- 
  LOGICAL, INTENT(out) :: errFound  !  TRUE indicates error detected

!------------------------------------------------------------------------------- 
! Array arguments with intent(out)
!------------------------------------------------------------------------------- 
  LOGICAL, INTENT(out) :: foundVar(:)            !  TRUE when a variable in the master
!                                                   list is found in the input list
!------------------------------------------------------------------------------- 
! Optional array arguments with intent(in)
!------------------------------------------------------------------------------- 
  INTEGER, INTENT(in), OPTIONAL :: varFlagIn(:)   !  flag for each variable
  INTEGER, INTENT(in), OPTIONAL :: varNzIn(:)     !  number of levels
  REAL, INTENT(in), OPTIONAL :: varConstIn(:)     !  constant values
  logical, INTENT(in), OPTIONAL :: varLog1In(:)   !  a logical switch
  logical, INTENT(in), OPTIONAL :: varLog2In(:)   !  a logical switch
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: varFileNameIn(:)  !  a character variable
  CHARACTER(LEN=*), INTENT(in), OPTIONAL :: varInterpIn(:)    !  a character variable
  CHARACTER(len=*), INTENT(in), OPTIONAL :: varNameSDFin(:)   !  name used for
!                                                  variable in a self-describing file

!------------------------------------------------------------------------------- 
! Optional array arguments with intent(out)
!------------------------------------------------------------------------------- 
  INTEGER, INTENT(out), OPTIONAL :: varFlag(:)  !  flag for each variable
  INTEGER, INTENT(out), OPTIONAL :: varNz(:)    !  number of levels
  REAL, INTENT(out), OPTIONAL :: varConst(:)    !  constant values
  logical, INTENT(out), OPTIONAL :: varLog1(:)  !  a logical switch
  logical, INTENT(out), OPTIONAL :: varLog2(:)  !  a logical switch
  CHARACTER(LEN=*), INTENT(out), OPTIONAL :: varFileName(:)  !  a character variable
  CHARACTER(LEN=*), INTENT(out), OPTIONAL :: varInterp(:)    !  a character variable
  CHARACTER(len=*), INTENT(out), OPTIONAL :: varNameSDF(:)   !  name used for
!                                                 variable in a self-describing file

!------------------------------------------------------------------------------- 
! Local scalar variables.
!------------------------------------------------------------------------------- 
  INTEGER :: ivar    !  loop counter
  INTEGER :: jvar    !  loop counter
  INTEGER :: nvar    !  size of input "master" list
  LOGICAL :: found   !  work

!------------------------------------------------------------------------------- 
! Initialise.
  errFound = .FALSE.
  nvar = SIZE( varNameList )
  foundVar(:) = .FALSE.  !  no variables in master list are in the input list
  IF ( PRESENT(varFlag) ) varFlag(:) = 0

!------------------------------------------------------------------ 

  DO ivar=1,nvarIn

    found = .FALSE.

!   Check if this variable can be found in the master list.
    DO jvar=1,nvar

      IF ( varNameIn(ivar) == varNameList(jvar) ) THEN
        found = .TRUE.

!       Load values for this variable.
        foundVar(jvar) = .TRUE.
        IF ( PRESENT(varFlagIn) ) THEN
          varFlag(jvar) = varFlagIn(ivar)
!         For an SDF, set varFlag to 1.
          IF ( sdfFile ) varFlag(jvar) = 1
        ENDIF
        IF ( PRESENT(varNameSDFin) ) varNameSDF(jvar) = varNameSDFin(ivar)
        IF ( PRESENT(varNzIn) ) varNz(jvar) = varNzIn(ivar)
        IF ( PRESENT(varConstIn) ) varConst(jvar) = varConstIn(ivar)
        IF ( PRESENT(varLog1In) ) varLog1(jvar) = varLog1In(ivar)
        IF ( PRESENT(varLog2In) ) varLog2(jvar) = varLog2In(ivar)
        IF ( PRESENT(varFileNameIn) ) varFileName(jvar) = varFileNameIn(ivar)
        IF ( PRESENT(varInterpIn) ) varInterp(jvar) = varInterpIn(ivar)

        EXIT  !  ivar variable is now done with
      ENDIF

    ENDDO

!   If a variable cannot be found, set error flag and print information.
    IF ( .NOT. found ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: varList: variable not recognised: ',TRIM(varNameIn(ivar))
    ENDIF

  ENDDO

! Extra information if error found.
  IF ( errFound ) THEN
    WRITE(*,*)'Error found by varList.'
    WRITE(*,*)'List of known variables:'
    DO ivar=1,nvar
      WRITE(*,*) TRIM( varNameList(ivar) )
    ENDDO
  ENDIF

  END SUBROUTINE varList

  END MODULE misc_utils
