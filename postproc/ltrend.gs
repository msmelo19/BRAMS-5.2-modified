***************************************************************************************
* $Id: ltrend.gs,v 1.18 2016/04/30 02:11:56 bguan Exp $
*
* Copyright (c) 2013-2015, Bin Guan
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are
* permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice, this list
*    of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice, this
*    list of conditions and the following disclaimer in the documentation and/or other
*    materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
* OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
* SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
* BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
***************************************************************************************
function ltrend(arg)
*
* Calculate linear trend over time (based on least-squares fitting).
*
rc=gsfallow('on')

input=subwrd(arg,1)
output=subwrd(arg,2)
slope=subwrd(arg,3)
rmse=subwrd(arg,4)
if(input=''); usage(); return; endif
if(output=''); output=input; endif

qdims(0,'mydim')

'set x '_.mydim.xs
'set y '_.mydim.ys
'set t '_.mydim.ts' '_.mydim.te
'tttmp=tloop(sum(1,t='_.mydim.ts',t+0))+'_.mydim.ts'-1'
'set x '_.mydim.xs' '_.mydim.xe
'set y '_.mydim.ys' '_.mydim.ye
'set t '_.mydim.ts
'trndslope=tregr(tttmp,'input',t='_.mydim.ts',t='_.mydim.te')'
'trndintercept=ave('input',t='_.mydim.ts',t='_.mydim.te')-trndslope*ave(maskout(tttmp,('input')-('input')+1),t='_.mydim.ts',t='_.mydim.te')'
'set t '_.mydim.ts' '_.mydim.te
'define trndrecon=trndslope*tttmp+trndintercept'
'define trnddiff=trndrecon-('input')'
if(output='DISPLAY' | output='display')
  'display trndrecon'
else
  'define 'output'=trndrecon'
endif
'set t '_.mydim.ts
if(slope!='')
  'define 'slope'=trndslope'
endif
if(rmse!='')
  'define 'rmse'=sqrt(ave(trnddiff*trnddiff,t='_.mydim.ts',t='_.mydim.te'))'
endif
'set t '_.mydim.ts' '_.mydim.te
'undefine tttmp'
'undefine trndslope'
'undefine trndintercept'
'undefine trndrecon'
'undefine trnddiff'

*
* Restore original dimension environment.
*
_.mydim.resetx
_.mydim.resety
_.mydim.resetz
_.mydim.resett

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Calculate linear trend over time (based on least-squares fitting).'
say ''
say '  USAGE: ltrend <input> [<output> [<slope> [<rmse>]]]'
say '    <input>: input field. Can be any GrADS expression.'
say '    <output>: output field, i.e., fitted trend line. Default=<input>.'
say '    <slope>: slope of fitted trend line, i.e., change of <input> per time step.'
say '    <rmse>: root mean square error.'
say ''
say '  DEPENDENCIES: qdims.gsf'
say ''
say '  Copyright (c) 2013-2015, Bin Guan.'
return
