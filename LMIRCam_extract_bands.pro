pro LMIRCam_extract_bands, phot=phot

path = '/home/amueller/work/LMIRCam/data/HD183324/Reduced/'

;=========================================================================

;PSF

if Keyword_set(phot) then begin

  filedx = file_search(path+'PSF_dx_*.fits', count=n)
  hdr = headfits(filedx[0],exten=0,/silent)
  radius = get_eso_keyword(hdr, 'NAXIS1')
  radius = (radius-1.)/2.

  tmpdx = dblarr(2*radius+1, 2*radius+1 ,n)
  for i=0,n-1 do tmpdx[*,*,i] = mrdfits(filedx[i],0,hdr,/silent)
  nframes = n_elements(tmpdx[0,0,*])
  medim = median(tmpdx, dim=3)
  sim = dblarr(nframes)
  for i=0,nframes-1 do sim[i] = stddev(tmpdx[*,*,i]-medim)
  resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
  spawn, 'rm '+path+'PSF_dx*fits'
  for i=0,n_elements(idxgood)-1 do writefits, path+'PSF_dx_'+strcompress(i+1,/rem)+'.fits', tmpdx[*,*,idxgood[i]], hdr
  psfdx = median(tmpdx[*,*,idxgood],dim=3,/even)
  writefits, path+'PSF_dx.fits', psfdx, hdr


  filesx = file_search(path+'PSF_sx_*.fits', count=n)
  tmpsx = dblarr(2*radius+1, 2*radius+1, n)
  for i=0,n-1 do tmpsx[*,*,i] = mrdfits(filesx[i],0,hdr,/silent)
  nframes = n_elements(tmpsx[0,0,*])
  medim = median(tmpsx, dim=3)
  sim = dblarr(nframes)
  for i=0,nframes-1 do sim[i] = stddev(tmpsx[*,*,i]-medim)
  resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
  spawn, 'rm '+path+'PSF_sx*fits'
  for i=0,n_elements(idxgood)-1 do writefits, path+'PSF_sx_'+strcompress(i+1,/rem)+'.fits', tmpsx[*,*,idxgood[i]], hdr
  psfsx = median(tmpsx[*,*,idxgood],dim=3,/even)
  writefits, path+'PSF_sx.fits', psfsx, hdr

  spawn, 'cp '+path+'PSF_dx.fits '+path+'PSF_sx.fits '+path+'../'
  
endif

;=========================================================================

;SX

file = file_search(path+'sx_cube_*reduced*', count=nfiles)
fileparal = file_search(path+'cube_*paral*')

dateim = strarr(nfiles)
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
pa = dblarr(nfiles)
ha = dblarr(nfiles)

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
  if (i eq 0) then begin

    nx = get_eso_keyword(hdr, 'NAXIS1')
    ny = nx
  
  endif
  
endfor

;-------------------------------------------------------------------------

;compute JD inside the exposure

for i=0,nfiles-1 do begin

  t1 = epoch[i]+'T'+tstart[i]
  t2 = epoch[i]+'T'+tend[i]

  jd1 = date_conv(t1, 'J')
  jd2 = date_conv(t2, 'J')

  jd[i] = mean([jd1,jd2])
  dateim[i] = t1

endfor

;-------------------------------------------------------------------------

idxu = uniq(roffx, sort(roffx))

final_sx = dblarr(nx,ny,n_elements(idxu))
final_pa = dblarr(n_elements(idxu))
final_ha = dblarr(n_elements(idxu))
final_jd = dblarr(n_elements(idxu))

for xx=0,n_elements(idxu)-1 do begin

  idx = where(roffx eq roffx[idxu[xx]])

  tmpim = dblarr(nx,ny,n_elements(idx))
  tmpha = dblarr(n_elements(idx))
  tmppa = dblarr(n_elements(idx))
  tmpjd = dblarr(n_elements(idx))
  
  for i=0,n_elements(idx)-1 do begin
  
    tmpim[*,*,i] = mrdfits(file[idx[i]],0,hdrsx,/silent)
    tmppa[i] = mrdfits(fileparal[idx[i]],0,/silent)
    tmpha[i] = mrdfits(fileparal[idx[i]],1,/silent)
    tmpjd[i] = jd[idx[i]]
    
  endfor
  
  if (n_elements(tmpim[0,0,*]) gt 1) then begin
  
    medim = median(tmpim, dim=3)
    sim = dblarr(n_elements(tmpim[0,0,*]))
    for i=0,n_elements(tmpim[0,0,*])-1 do sim[i] = stddev(tmpim[*,*,i]-medim)
    resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
    final_sx[*,*,xx] = median(tmpim[*,*,idxgood], dim=3, /even)
    final_pa[xx] = mean(tmppa[idxgood])
    final_ha[xx] = mean(tmpha[idxgood])
    final_jd[xx] = mean(tmpjd[idxgood])

  endif else begin
  
    final_sx[*,*,xx] = tmpim
    final_pa[xx] = tmppa
    final_ha[xx] = tmpha
    final_jd[xx] = tmpjd
  
  endelse
    
endfor

idxs = sort(final_jd)
final_sx = final_sx[*,*,idxs]
final_pa_sx = final_pa[idxs]
final_ha_sx = final_ha[idxs]

;=========================================================================

;DX

file = file_search(path+'dx_cube_*reduced*', count=nfiles)
fileparal = file_search(path+'cube_*paral*')

dateim = strarr(nfiles)
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
pa = dblarr(nfiles)
ha = dblarr(nfiles)

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
  if (i eq 0) then begin

    nx = get_eso_keyword(hdr, 'NAXIS1')
    ny = nx
  
  endif
  
endfor

;-------------------------------------------------------------------------

;compute JD inside the exposure

for i=0,nfiles-1 do begin

  t1 = epoch[i]+'T'+tstart[i]
  t2 = epoch[i]+'T'+tend[i]

  jd1 = date_conv(t1, 'J')
  jd2 = date_conv(t2, 'J')

  jd[i] = mean([jd1,jd2])
  dateim[i] = t1

endfor

;-------------------------------------------------------------------------

idxu = uniq(roffx, sort(roffx))

final_dx = dblarr(nx,ny,n_elements(idxu))
final_pa = dblarr(n_elements(idxu))
final_ha = dblarr(n_elements(idxu))
final_jd = dblarr(n_elements(idxu))

for xx=0,n_elements(idxu)-1 do begin

  idx = where(roffx eq roffx[idxu[xx]])

  tmpim = dblarr(nx,ny,n_elements(idx))
  tmpha = dblarr(n_elements(idx))
  tmppa = dblarr(n_elements(idx))
  tmpjd = dblarr(n_elements(idx))
  
  for i=0,n_elements(idx)-1 do begin
  
    tmpim[*,*,i] = mrdfits(file[idx[i]],0,hdrdx,/silent)
    tmppa[i] = mrdfits(fileparal[idx[i]],0,/silent)
    tmpha[i] = mrdfits(fileparal[idx[i]],1,/silent)
    tmpjd[i] = jd[idx[i]]
    
  endfor
  
  if (n_elements(tmpim[0,0,*]) gt 1) then begin
  
    medim = median(tmpim, dim=3)
    sim = dblarr(n_elements(tmpim[0,0,*]))
    for i=0,n_elements(tmpim[0,0,*])-1 do sim[i] = stddev(tmpim[*,*,i]-medim)
    resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgood, badvec=idxbad
    final_dx[*,*,xx] = median(tmpim[*,*,idxgood], dim=3, /even)
    final_pa[xx] = mean(tmppa[idxgood])
    final_ha[xx] = mean(tmpha[idxgood])
    final_jd[xx] = mean(tmpjd[idxgood])

  endif else begin
  
    final_dx[*,*,xx] = tmpim
    final_pa[xx] = tmppa
    final_ha[xx] = tmpha
    final_jd[xx] = tmpjd
  
  endelse

endfor

idxs = sort(final_jd)
final_dx = final_dx[*,*,idxs]
final_pa_dx = final_pa[idxs]
final_ha_dx = final_ha[idxs]

;=========================================================================

;final QC
medim = median(final_sx, dim=3)
sim = dblarr(n_elements(final_sx[0,0,*]))
for i=0,n_elements(final_sx[0,0,*])-1 do sim[i] = stddev(final_sx[*,*,i]-medim)
resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgoodsx, badvec=idxbadsx

medim = median(final_dx, dim=3)
sim = dblarr(n_elements(final_dx[0,0,*]))
for i=0,n_elements(final_dx[0,0,*])-1 do sim[i] = stddev(final_dx[*,*,i]-medim)
resistant_mean, sim, 3, t1, t2, nbad, /double, goodvec=idxgooddx, badvec=idxbaddx

writefits, path+'img_sx_dc.fits', final_sx[*,*,idxgoodsx], hdrsx
writefits, path+'vec_sx_paral.fits', final_pa_sx[idxgoodsx]
writefits, path+'vec_sx_paral.fits', final_ha_sx[idxgoodsx], /append

writefits, path+'img_dx_dc.fits', final_dx[*,*,idxgooddx], hdrdx
writefits, path+'vec_dx_paral.fits', final_pa_dx[idxgooddx]
writefits, path+'vec_dx_paral.fits', final_ha_dx[idxgooddx], /append

spawn, 'cp '+path+'img*fits '+path+'vec*fits '+path+'../'

end
