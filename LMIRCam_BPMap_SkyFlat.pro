;flat made out of science frames themselve

pro LMIRCam_BPMap_SkyFlat

;x0,x1,y0,y1
;box = [382,2047,59,1279]
box = [0,2047,0,1279]
;assuming /xxx/yyy/<star>/RAW
path = ''
;read, 'Enter full path to raw data: ', path
path = '/home/amueller/work/LMIRCam/data/HD183324/FF/'
path = path+'/'

resdir = path+'../Reduced/'
spawn, 'mkdir -p '+resdir

file = file_search(path+'*fits', count=nfiles)

;-------------------------------------------------------------------------

tdit = dblarr(nfiles)

for i=0,nfiles-1 do begin

  hdr = headfits(file[i], exten=0, /silent)
  tdit[i] = double(get_eso_keyword(hdr, 'ITIME'))
  if (i eq 0) then begin

    nx = get_eso_keyword(hdr, 'NAXIS1')
    ny = get_eso_keyword(hdr, 'NAXIS2')

  endif

endfor

idxs = sort(tdit)
file = file[idxs]
tdit = tdit[idxs]

tmp = dblarr(nx, ny, nfiles)
for i=0,nfiles-1 do tmp[*,*,i] = mrdfits(file[i],0,/silent,/dscale)

tmp = tmp[box[0]:box[1],box[2]:box[3],*]
nx = n_elements(tmp[*,0,0])
ny = n_elements(tmp[0,*,0])

idxu = uniq(tdit)
ff = dblarr(nx,ny,n_elements(idxu))
dit = tdit[idxu]
sdev = dblarr(n_elements(idxu))
for i=0,n_elements(idxu)-1 do begin

  idx = where(tdit eq dit[i])
  ff[*,*,i] = mean(tmp[*,*,idx], dim=3)

endfor

tmp = 0

;-------------------------------------------------------------------------

slope = dblarr(nx,ny)

for xx=0,nx-1 do begin

  for yy=0,ny-1 do begin

    sixlin, dit, ff[xx,yy,*], a, siga, b, sigb
    slope[xx,yy] = b[0]

  endfor

  proceeding_text,loop=(nx), i=xx, prompt='> x direction   '+string(xx+1,form='(I4)')

endfor

;-------------------------------------------------------------------------

; /home/amueller/work/LMIRCam/data/HD183324/FF
bp = intarr(nx,ny)
bp[*,*] = 1
rm = robust_mean(slope, 5., sigma, numrej, goodind=good)

idxgood = array_indices(bp, good)
for i=0L,n_elements(idxgood[0,*])-1 do bp[idxgood[0,i], idxgood[1,i]] = 0

writefits, resdir+'BadPixelMap.fits', bp

;-------------------------------------------------------------------------


;normalisation

;median slope
idx = where(bp eq 0.d0)
msl = median(slope[idx], /even)

ff_final = slope/msl

writefits, resdir+'FlatField.fits', ff_final

end
