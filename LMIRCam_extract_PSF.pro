pro LMIRCam_extract_PSF

;x0,x1,y0,y1
; box = [382,2047,59,1279]
box = [0,2047,0,1279]

; path = ''
; read, 'Enter full path to raw data: ', path
; path = path+'/'
path = '/home/amueller/work/LMIRCam/data/HD183324/Flux/'

plsc = 10.707d-3
diam = 8.4
lambda = 3.7d-6	;L'
fwhm = lambda/diam*206265.d0/plsc

resdir = path+'../Reduced/'
spawn, 'mkdir -p '+resdir

qintsky = 'n'

;-------------------------------------------------------------------------

bp = mrdfits(resdir+'BadPixelMap.fits', 0, /silent)
ff = mrdfits(resdir+'FlatField.fits', 0, /silent)

nx = n_elements(bp[*,0])
ny = n_elements(bp[0,*])

;-------------------------------------------------------------------------

skf = file_search(resdir+'PSF_sky*background*fits', count=nskf)

pos1 = strpos(skf, 'PSF_sky_')+8
step = strmid(skf, pos1[0], 1)
step = step[uniq(step)]

pos2 = strpos(skf, '_', /reverse_search)
pos3 = strpos(skf, '.fits', /reverse_search)
skyreloffs = strarr(nskf)
for i=0,nskf-1 do skyreloffs[i] = strmid(skf[i], pos2[i]+1, pos3[i]-pos2[i]-1)
skyureloffs = skyreloffs[uniq(skyreloffs, sort(skyreloffs))]

;-------------------------------------------------------------------------

for qu=0,n_elements(step)-1 do begin
;for qu=0,0 do begin

  ;-----------------------------------------------------------------------
  ;find the right frames
  file = file_search(path+'*fits', count=nfiles)

  jd = dblarr(nfiles)
  roffx = dblarr(nfiles)
  roffy = dblarr(nfiles)
  loffx = dblarr(nfiles)
  loffy = dblarr(nfiles)
  dit = dblarr(nfiles)
  epoch = strarr(nfiles)
  tstart = strarr(nfiles)
  tend = strarr(nfiles)
  dateim = strarr(nfiles)

  for i=0,nfiles-1 do begin

    hdr = headfits(file[i], exten=0, /silent)
    roffx[i] = get_eso_keyword(hdr, 'ROFFSETX')
    roffy[i] = get_eso_keyword(hdr, 'ROFFSETY')
    loffx[i] = get_eso_keyword(hdr, 'LOFFSETX')
    loffy[i] = get_eso_keyword(hdr, 'LOFFSETY')
    epoch[i] = strcompress(get_eso_keyword(hdr, 'DATE-OBS'))
    tstart[i] = strcompress(get_eso_keyword(hdr, 'TIME-OBS'))
    tend[i] = strcompress(get_eso_keyword(hdr, 'TIME-END'))
    dit[i] = double(get_eso_keyword(hdr, 'ITIME'))

  endfor

  ;compute JD inside the exposure

  for i=0,nfiles-1 do begin

    t1 = epoch[i]+'T'+tstart[i]
    t2 = epoch[i]+'T'+tend[i]

    jd1 = date_conv(t1, 'J')
    jd2 = date_conv(t2, 'J')

    jd[i] = mean([jd1,jd2])
    dateim[i] = t1

  endfor

  signx = sign(roffx)
  signy = sign(roffy)
  reloffs = strcompress(signx,/rem)+strcompress(signy,/rem)
  ;ureloffs = reloffs[uniq(reloffs, sort(reloffs))]

  idx = where(reloffs ne skyureloffs[qu])	;'ne' otherwise sky has the some offset
  rawf = file[idx]
  date = dateim[idx]
  jd = jd[idx]

  ;-----------------------------------------------------------------------


  ;idx = where(flag eq qu+1)
  ;rawf = path+rfname[idx]

  ;tmp = headfits(rawf[0],exten=0,/silent)
  ;nx = get_eso_keyword(tmp, 'NAXIS1')
  ;ny = get_eso_keyword(tmp, 'NAXIS2')


  ;sky background
  skf = file_search(path+'../Reduced/'+'PSF_sky_'+strcompress(qu+1,/rem)+'_background_*_'+skyureloffs[qu]+'.fits', count=nsk)

  ;extract date from sky files to find the closest one to science later
  jdsk = dblarr(n_elements(skf))
  for i=0,n_elements(skf)-1 do begin

    head = headfits(skf[i], exten=0, /silent)
    jdsk[i] = double(get_eso_keyword(head, 'JD'))

  endfor

  ;read in all sky frames for later interpolation
  hdrtmp = headfits(skf[0], exten=0, /silent)
  sknx = float(get_eso_keyword(hdrtmp, 'NAXIS1'))
  skny = float(get_eso_keyword(hdrtmp, 'NAXIS2'))
  sktmp = dblarr(sknx, skny, n_elements(skf))
  nsky = n_elements(skf)
  for j=0,nsky-1 do sktmp[*,*,j] = mrdfits(skf[j], 0, /silent)

  dimx = sknx
  dimy = skny

  ;==================================================================================================

  ;read in cubes, do cosmetics, average

  nfiles = n_elements(rawf)
  ;jd = dblarr(nfiles)
  hafile = dblarr(nfiles)
  used_frames = intarr(2,nfiles)
  outname = strarr(nfiles)

  flagcut = intarr(nfiles)

  ;==================================================================================================


  for i=0,nfiles-1 do begin
  ;for i=55,55 do begin

    print, ''
    print, 'Reducing Cube '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
    print, 'Reading File'

    cube = mrdfits(rawf[i], 0, hdrraw, /silent, /dscale)
    cube = cube[box[0]:box[1],box[2]:box[3]]

    szcube = size(cube)
    nframes = n_elements(cube[0,0,*])

    ;----------------------------------------------------------------------------------

    ;output name
    pos1 = strpos(rawf[i], '/', /reverse_search)
    pos2 = strpos(rawf[i], '.fits', /reverse_search)
    outname[i] = strmid(rawf[i], pos1+1, pos2-pos1-1)

    ;----------------------------------------------------------------------------------

    ;find/create corresponding sky frame

    if (qintsky eq 'n') then begin

      idx = closest(jdsk, jd[i])
      sk = mrdfits(skf[idx],0,/silent)

    endif

    if (qintsky eq 'y') then begin

      if (n_elements(jdsk) gt 1) then begin

	idx = closest(jdsk, jd[i])

	if (jd[i] lt jdsk[0]) then sk = mrdfits(skf[idx],0,/silent)	;first obs with sky afterwards, no interpolation needed

	if (jd[i] gt jdsk[n_elements(jdsk)-1]) then sk = mrdfits(skf[idx],0,/silent)	;last observation is target and no sky, again no interpolation
	  
	if (jd[i] gt jdsk[0] and jd[i] lt jdsk[n_elements(jdsk)-1]) then begin

	  sk = dblarr(nx,ny)

	  for xx=0,sknx-1 do begin
	    for yy=0,skny-1 do begin
	      sk[xx,yy] = interpol(sktmp[xx,yy,*], jdsk, jd[i])
	    endfor
	  endfor

	endif

	;sk = sk[*,0:nx-1]

      endif else begin

	sk = mrdfits(skf[0],0,/silent)
	;sk = sk[*,0:nx-1]

      endelse

    endif

    ;----------------------------------------------------------------------------------

    im = fltarr(nx,ny,nframes)
    corim = fltarr(nx,ny,nframes)	;cosmetically corrected images
    sdev = dblarr(nframes)
    flux = dblarr(nframes)
    mflux = dblarr(nframes)

    ;----------------------------------------------------------------------------------

    for j=0,nframes-1 do begin

      ;cosmetics

      if (qintsky eq 'y' or qintsky eq 'n') then begin

	im[*,*,j] = (cube[*,*,j]-sk)/ff
	tmp = im[*,*,j]
	fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
	corim[*,*,j] = outim

      endif

      ;----------------------------------------------------------------------------------


    endfor	;nframes

    ;----------------------------------------------------------------------------------


    ;find stars
    synthg = gaussian2d(80,80,40,40,fwhm,fwhm, 0)
    starfinder, corim, synthg, [max(corim)/4., max(corim)/4.], 0.7, xpos, ypos, fluxes, sx, sy, sf, correlation, /silent

    if (n_elements(xpos) gt 2.) then begin

      print, ''
      print, '********************************'
      print, 'More than 2 stars detected.'
      print, '********************************'
      stop

    endif

    ;if only one star is detected insert dummy values because AO might have stopped
    if (n_elements(xpos) eq 1.) then begin

      if (ypos lt 650.) then begin

	xpos = [xpos,xpos]
	ypos = [ypos,900.]

      endif else begin

	xpos = [xpos,xpos]
	ypos = [ypos,350.]

      endelse

    endif

    ;cut images for dx (up) and sx (down)

    radius = 22.

    ;SX
    idx = where(ypos lt 650.)
    xc = floor(xpos[idx])-xpos[idx]
    yc = floor(ypos[idx])-ypos[idx]
    stmp = fftshift(corim, xc, yc)
    s1x = floor(xpos[idx])-radius
    s2x = floor(xpos[idx])+radius
    s1y = floor(ypos[idx])-radius
    s2y = floor(ypos[idx])+radius
    flag_sx = 0

    if (s1x ge 0. and s2x lt (size(stmp,/dim))[0] and s1y ge 0. and s2y lt (size(stmp,/dim))[1]) then begin
    
      cutim_sx = stmp[s1x:s2x,s1y:s2y]
      flag_sx = 1
    
    endif else begin
    
      if (s1x lt 0.) then s1x = 0
      if (s1y lt 0.) then s1y = 0
      if (s2x ge (size(stmp,/dim))[0]) then s2x = (size(stmp,/dim))[0]-1.
      if (s2y ge (size(stmp,/dim))[1]) then s2y = (size(stmp,/dim))[1]-1.
    
    endelse
    

    ;DX
    idx = where(ypos gt 650.)
    xc = floor(xpos[idx])-xpos[idx]
    yc = floor(ypos[idx])-ypos[idx]
    stmp = fftshift(corim, xc, yc)
    s1x = floor(xpos[idx])-radius
    s2x = floor(xpos[idx])+radius
    s1y = floor(ypos[idx])-radius
    s2y = floor(ypos[idx])+radius
    flag_dx = 0

    if (s1x ge 0. and s2x lt (size(stmp,/dim))[0] and s1y ge 0. and s2y lt (size(stmp,/dim))[1]) then begin
    
      cutim_dx = stmp[s1x:s2x,s1y:s2y]
      flag_dx = 1
    
    endif else begin
    
      if (s1x lt 0.) then s1x = 0
      if (s1y lt 0.) then s1y = 0
      if (s2x ge (size(stmp,/dim))[0]) then s2x = (size(stmp,/dim))[0]-1.
      if (s2y ge (size(stmp,/dim))[1]) then s2y = (size(stmp,/dim))[1]-1.
    
    endelse

    ;output
    if (flag_dx eq 1) then writefits, path+'../Reduced/'+'PSF_dx_'+strcompress(i+1,/rem)+'.fits', cutim_dx, hdrraw
    if (flag_sx eq 1) then writefits, path+'../Reduced/'+'PSF_sx_'+strcompress(i+1,/rem)+'.fits', cutim_sx, hdrraw
    
    ;proceeding_text,loop=nfiles, i=i, prompt='> Processing Frame   '+string(i+1,form='(I4)')

  endfor	;nfiles

endfor


stop
end
