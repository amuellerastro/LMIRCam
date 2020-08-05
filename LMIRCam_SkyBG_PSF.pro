pro LMIRCam_SkyBG_PSF

;x0,x1,y0,y1
;box = [382,2047,59,1279]
box = [0,2047,0,1279]

;assuming /xxx/yyy/<star>/RAW
path = ''
;read, 'Enter full path to raw data: ', path
path = '/home/amueller/work/LMIRCam/data/HD183324/Flux'
path = path+'/'

resdir = path+'../Reduced/'
spawn, 'mkdir -p '+resdir

rawf = file_search(path+'*fits', count=nfiles)

;-------------------------------------------------------------------------
; https://zero.as.arizona.edu/wiki/pages/H0v8N3Q3Y/Pixel_Scale_and_Orientation.html
plsc = 10.707d-3
diam = 8.4
lambda = 3.7d-6	;L'
fwhm = lambda/diam*206265.d0/plsc

;-------------------------------------------------------------------------

bp = mrdfits(resdir+'BadPixelMap.fits', 0, /silent)

pos = strpos(rawf, '/', /reverse_search)
rfname = strmid(rawf, pos[0]+1, strlen(rawf[0])-pos[0]-1)

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

for i=0,nfiles-1 do begin

  hdr = headfits(rawf[i], exten=0, /silent)
  roffx[i] = get_eso_keyword(hdr, 'ROFFSETX')
  roffy[i] = get_eso_keyword(hdr, 'ROFFSETY')
  loffx[i] = get_eso_keyword(hdr, 'LOFFSETX')
  loffy[i] = get_eso_keyword(hdr, 'LOFFSETY')
  epoch[i] = strcompress(get_eso_keyword(hdr, 'DATE-OBS'))
  tstart[i] = strcompress(get_eso_keyword(hdr, 'TIME-OBS'))
  tend[i] = strcompress(get_eso_keyword(hdr, 'TIME-END'))
  dit[i] = double(get_eso_keyword(hdr, 'ITIME'))

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

signx = sign(roffx)
signy = sign(roffy)
reloffs = strcompress(signx,/rem)+strcompress(signy,/rem)
;reloffs = strcompress(roffx,/rem)+strcompress(roffy,/rem)
ureloffs = reloffs[uniq(reloffs, sort(reloffs))]

nu = n_elements(ureloffs)

for i=0,nu-1 do begin

  idx = where(reloffs eq ureloffs[i])

  d = ts_diff(idx,1)
  idx1 = where(d ne -1.)

  for j=0,n_elements(idx1)-1 do begin

    if (j eq 0) then begin

      tmpim = dblarr(n_elements(bp[*,0]),n_elements(bp[0,*]),idx1[0]+1)
      tmpjd = dblarr(idx1[0]+1)
      tmpdate = strarr(idx1[0]+1)

    endif

    if (j ne 0) then begin

      tmpim = dblarr(n_elements(bp[*,0]),n_elements(bp[0,*]),idx1[j]-idx1[j-1])
      tmpjd = dblarr(idx1[j]-idx1[j-1])
      tmpdate = strarr(idx1[j]-idx1[j-1])

    endif

    for k=0,n_elements(tmpim[0,0,*])-1 do begin

      if (j eq 0) then begin

	if (k eq round(n_elements(tmpim[0,0,*])/2.)) then tmp = mrdfits(rawf[idx[k]],0,hdr,/silent,/dscale) else tmp = mrdfits(rawf[idx[k]],0,/silent,/dscale)

	tmpjd[k] = jd[idx[k]]
	tmpdate[k] = dateim[idx[k]]

      endif

      if (j ne 0) then begin

	if (k eq round(n_elements(tmpim[0,0,*])/2.)) then tmp = mrdfits(rawf[idx[idx1[j-1]+1+k]],0,hdr,/silent,/dscale) else tmp = mrdfits(rawf[idx[idx1[j-1]+1+k]],0,/silent,/dscale)
	tmpjd[k] = jd[idx[idx1[j-1]+1+k]]
	tmpdate[k] = dateim[idx[idx1[j-1]+1+k]]

      endif

      tmp = tmp[box[0]:box[1],box[2]:box[3]]
      tmpim[*,*,k] = tmp

    endfor

    sky = mean(tmpim, dim=3)
    tmp = sky
    fixpix_mod, tmp, bp, outim, npix=24, /weight, /silent
    sky = outim
    sxaddpar, hdr, 'JD', mean(tmpjd), format='(f18.10)'

    skyname = path+'../Reduced/'+'PSF_sky_'+strcompress(i+1,/rem)+'_background_'+tmpdate[round(n_elements(tmpim[0,0,*])/2.)]+'_'+ureloffs[i]+'.fits'
    writefits, skyname, outim, hdr

  endfor

  ;proceeding_text,loop=(i), i=nu, prompt='> Sky   '+string(i+1,form='(I4)')

endfor


stop
end
