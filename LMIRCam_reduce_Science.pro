@eq2hor_MOD.pro

function get_parangle, jd, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres, hafile

  caldat, jd, mo, day, yr, hh, mm, ss

  ;get date of the form yyyy-mm-ddThh:mm:ss
  date = dblarr(6)
  date[0] = yr
  date[1] = mo
  date[2] = day
  date[3] = hh
  date[4] = mm
  date[5] = ss

  ;get year of observations for new RA, DEC, i.e. epoch
  hour = hh+mm/60.d0+ss/3600.d0
  if ((double(yr) mod 4.) eq 0.) then date2 = yr+mo/12.d0+day/366.d0+hour/8766.0d0 $
    else date2 = yr+mo/12.d0+day/365.d0+hour/8766.0d0	;takes leap year into account

  ;compute current cordinates and parallactic angle at time of observation
  dt = date2-epoch
  eq2hor_MOD, ra, dec, pma, pmd, radvel, plx, dt, jd, alt, az, hatmp, lat=lat, lon=lon, altitude=altitude, pres=pres, temp=temp, outra=outra, outdec=outdec;, /verbose
  ha = hatmp
  ha = ha/15.

; outdec = -41.39791667
  parang = parangle(ha, outdec, lat)

  if (outdec gt lat and parang lt 0.) then parang = parang+360.d0

  return, {parang:parang, ha:ha}
;   return, {parang:parang}

end

pro LMIRCam_reduce_Science

;x0,x1,y0,y1
; box = [382,2047,59,1279]
box = [0,2047,0,1279]

;assuming /xxx/yyy/<star>/RAW
; path = ''
; read, 'Enter full path to raw data: ', path
; path = path+'/'
path = '/home/amueller/work/LMIRCam/data/HD183324/RAW/'
filestar = '/home/amueller/work/IDLlibs/AO/TargetProperties/Targets/HD183324.sav'
restore, filestar, /verbose
ra_st = ra
dec_st = dec

plsc = 10.707d-3
diam = 8.4
lambda = 3.7d-6	;L'
fwhm = lambda/diam*206265.d0/plsc

resdir = path+'../Reduced/'
spawn, 'mkdir -p '+resdir

; qintsky = ''
; read, 'SkyPCA (p) or Interpolate Sky frames (y/n): ', qintsky
qintsky = 'n'
;PCA does not work yets

;-------------------------------------------------------------------------

bp = mrdfits(resdir+'BadPixelMap.fits', 0, /silent)
ff = mrdfits(resdir+'FlatField.fits', 0, /silent)

nx = n_elements(bp[*,0])
ny = n_elements(bp[0,*])

;-------------------------------------------------------------------------

skf = file_search(resdir+'sky*background*fits', count=nskf)

pos1 = strpos(skf, 'sky_')+4
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
  skf = file_search(path+'../Reduced/'+'sky_'+strcompress(qu+1,/rem)+'_background_*_'+skyureloffs[qu]+'.fits', count=nsk)

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
  outname = strarr(nfiles)

  ;==================================================================================================

  ;load all sky frames and compute PCs and Eigenvalues for later

  if (qintsky eq 'p') then begin

    truncate_pca = nsky
    obj = sktmp
    nobj = (size(obj))[3]

    for i=0,nobj-1 do begin

      dum = 0
      dum = reform(obj[0:dimx-1,0:dimy-1,i])
      dum = dum-mean(dum)
      obj[0:dimx-1,0:dimy-1,i] = dum

    endfor

    ;PCA
    data = transpose(reform(obj,dimx*dimy,nobj))	;[0:dim-1,0:dim-1,*]
    covMatrix = matrix_multiply(data, data, /btranspose)

    eigenval = la_eigenql(covMatrix, EIGENVECTORS=eigenvect, range=[nobj-truncate_pca,nobj-1], /DOUBLE)
    eigenval = reverse(eigenval)

    eigenvect = reverse(eigenvect,2)
    pc_orig = matrix_multiply(eigenvect,data,/atranspose)
    pc = pc_orig
    for k=0,nsky-1 do pc[k,*] = pc_orig[k,*]/(eigenval[k])

  endif

  ;==================================================================================================

  for i=0,nfiles-1 do begin
  ;for i=55,55 do begin

    print, ''
    print, 'Reducing Cube '+strcompress(i+1,/rem)+' / '+strcompress(nfiles,/rem)
    print, 'Reading File'

    cube = mrdfits(rawf[i], 0, hdrraw, /silent, /dscale)
    cube = cube[box[0]:box[1],box[2]:box[3]]


    if (i eq 0) then hdrsciref = hdrraw

    szcube = size(cube)
;     if (szcube[1] lt szcube[2]) then cube = cube[*,0:szcube[1]-1,*]	;assumption: X additional rows in y-dir

;     if (n_elements(szcube) gt 5) then cube = cube[*,*,1:n_elements(cube[0,0,*])-2]	;last frame is already averaged, remove it, 1st frame looks bad as well

    nframes = n_elements(cube[0,0,*])
  ;   cube = cube;+32768.

    ;----------------------------------------------------------------------------------

    ;Observatory parameters
    lat = ten(+32.d0, 42.d0, 00.0d0)	;
    lon = ten(-109.d0, 50.d0, 59.99d0)	;
    altitude = 3170.d0

    ;----------------------------------------------------------------------------------

    ;output name
    pos1 = strpos(rawf[i], '/', /reverse_search)
    pos2 = strpos(rawf[i], '.fits', /reverse_search)
    outname[i] = strmid(rawf[i], pos1+1, pos2-pos1-1)

    lbt_lst = get_eso_keyword(hdrraw, 'LBT_LST')
    lbt_ra = get_eso_keyword(hdrraw, 'LBT_RA')

    t1 = strmid(lbt_lst, 0, strpos(lbt_lst,':'))
    t2 = strmid(lbt_lst, strpos(lbt_lst,':')+1, 2)
    t3 = strmid(lbt_lst, strpos(lbt_lst,':',/reverse_search)+1, 6)
    lbt_lst = ten(t1, t2, t3)

    t1 = strmid(lbt_ra, 0, strpos(lbt_ra,':'))
    t2 = strmid(lbt_ra, strpos(lbt_ra,':')+1, 2)
    t3 = strmid(lbt_ra, strpos(lbt_ra,':',/reverse_search)+1, 6)
    lbt_ra = ten(t1, t2, t3)

    hafile[i] = lbt_lst-lbt_ra

    ;pres = double(get_eso_keyword(hdrraw, 'HIERARCH ESO TEL AMBI PRES START'))
    ;pres = pres*1.33322	;mmHg to mbar
    temp = double(get_eso_keyword(hdrraw, 'LBTTEMP'))+273.15d0
    epoch = 2000.

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
    ;     im[*,*,j] = (cube[*,*,j]-(sk*(median(cube[*,*,j])/median(sk))))/ff	;subtract sky
	tmp = im[*,*,j]
	fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
	corim[*,*,j] = outim

      endif

      if (qintsky eq 'p') then begin

	tim = cube[*,*,j]
	fixpix_mod, tim, bp, outim, npix=24, /weight, /silent
	outim = outim-mean(outim)

	data2 = transpose(reform(outim,dimx*dimy))
	s1 = matrix_multiply(pc_orig, data2, /btranspose)
	sk = reform(matrix_multiply(s1, pc), dimx, dimy)

	corim[*,*,j] = (outim-sk)/ff

      endif

      ;----------------------------------------------------------------------------------

;       ;remove additional horizontal additive pattern by subtracting the median of each row (NACO manual, p62 or Sec.8.6: Cube mode)
;       ;we have 4 quadrants, hence define windows for each quadrant
; 
;       if (nx gt 500) then ws = 40 else ws = 20
;       ll = median(corim[0:ws-1,*,j], dimension=1)	;left det., left window
;       lr = median(corim[(nx/2)-ws:(nx/2)-1,*,j], dimension=1)	;left det., right window
;       rl = median(corim[nx/2:nx/2+ws-1,*,j], dimension=1)	;right det., left window
;       rr = median(corim[nx-1-ws:*,*,j], dimension=1)	;right det., right window
; 
;       for k=0,ny-1 do begin
; 
; 	corim[0:nx/2-1,k,j] = corim[0:nx/2-1,k,j]-mean([ll[k],lr[k]])
; 	corim[nx/2:*,k,j] = corim[nx/2:*,k,j]-mean([rl[k],rr[k]])
; 
;       endfor


      ;----------------------------------------------------------------------------------


  ;     ;mask for quality control if jitter is used
  ;     if (qagpm eq 'n') then begin
  ; 
  ;       txc = nx/2.+offx	;estimated positions
  ;       tyc = ny/2.+offy
  ; 
  ;       mask_t = shift(dist(nx), txc, tyc)
  ;       mask = mask_t ge 1.*fwhm and mask_t le 4.*fwhm
  ;       mask = mask[*]	;convert to 1D vector
  ;       idxmask = where(mask eq 1)
  ; 
  ;     endif

    endfor	;nframes

;     dum = max(corim, idxmax)
;     idx0 = array_indices(corim, idxmax)
;     ;estimated OBJECT positions
;     txc = idx0[0]
;     tyc = idx0[1]

    ;----------------------------------------------------------------------------------

    ;compute parallactic angle for exposure start and end of a cube

    ra = [double(strmid(ra_st,0,2)),double(strmid(ra_st,2,2)),double(strmid(ra_st,4,6))]
    sign1 = strmid(dec_st,0,1)
    if (sign1 eq '-') then sign1 = -1.d0 else sign1 = 1.d0
    dec = [sign1*double(strmid(dec_st,1,2)), double(strmid(dec_st,3,2)), double(strmid(dec_st,5,6))]
    ra = ten(ra)*15.
    dec = ten(dec)

    if (finite(radvel[0] ne 1)) then radvel = 0.
    if (finite(plx[0] ne 1)) then plx = 1./100.

;     if (n_elements(size(corim)) gt 5) then begin
; 
;       ;+1 because I remove the first frame anyway
;       st = get_parangle(jd[i]+(dit*(idxgood[0]+1.))/86400.d0, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
;       parang_start = st.parang
; 
;       st = get_parangle(jd[i]+(dit*(idxgood[n_elements(idxgood)-1]+1.))/86400.d0, epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres)
;       parang_end = st.parang
; 
;       parang = (parang_start+parang_end)/2.
; 
;     endif else begin

      st = get_parangle(jd[i], epoch, ra, dec, pma, pmd, radvel, plx, lat, lon, altitude, temp, pres, hafile[i])
      parang = st.parang


;     endelse

    ha = st.ha
    if (ha gt 12.) then ha = ha-24.
    if (ha lt -12.) then ha = ha+24.

    PA_onsky = parang

    ;----------------------------------------------------------------------------------

    ;find stars
    synthg = gaussian2d(80,80,40,40,fwhm,fwhm, 0)
    starfinder, corim, synthg, [max(corim)/2., max(corim)/2.], 0.7, xpos, ypos, fluxes, sx, sy, sf, correlation, /silent

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

    radius = 200.

    ;SX
    idx = where(ypos lt 650.)
    xc = floor(xpos[idx])-xpos[idx]
    yc = floor(ypos[idx])-ypos[idx]
    stmp = fftshift(corim, xc, yc)
    s1x = floor(xpos[idx])-radius
    s2x = floor(xpos[idx])+radius-1
    s1y = floor(ypos[idx])-radius
    s2y = floor(ypos[idx])+radius-1
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
    s2x = floor(xpos[idx])+radius-1
    s1y = floor(ypos[idx])-radius
    s2y = floor(ypos[idx])+radius-1
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
    if (flag_dx eq 1) then writefits, path+'../Reduced/'+'dx_cube_'+outname[i]+'_reduced.fits', cutim_dx, hdrraw
    if (flag_sx eq 1) then writefits, path+'../Reduced/'+'sx_cube_'+outname[i]+'_reduced.fits', cutim_sx, hdrraw
    if (flag_dx eq 1 and flag_sx eq 1) then begin

      writefits, path+'../Reduced/'+'cube_'+outname[i]+'_paral.fits', [PA_onsky]
      writefits, path+'../Reduced/'+'cube_'+outname[i]+'_paral.fits', [ha], /append

    endif
    
    ;proceeding_text,loop=nfiles, i=i, prompt='> Processing Frame   '+string(i+1,form='(I4)')

  endfor	;nfiles

endfor


stop
end
