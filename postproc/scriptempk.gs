* gera imagem da temperatura
'clear'
'set mpdset brmap_hires'
'open METEO-ONLY-A-2015-08-31-000000-g1.ctl'
*'paletatemp.gs'
'BTR01.gs'
'set gxout shaded'
'd tempk-273'
'set gxout contour'
'set cint 8'
'd tempk-273'
'draw title Temperatura'
*'draw string '
'cbarn.gs'
'printim medium_size.png white'
*'gxprint medium_size.png white'
'set rbcols'


