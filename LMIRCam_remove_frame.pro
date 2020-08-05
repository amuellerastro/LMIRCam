@readcol.pro
@remchar.pro
@gettok.pro
@strsplit.pro
@mrdfits.pro
@fxposit.pro
@mrd_hread.pro
@fxpar.pro
@valid_num.pro
@mrd_skip.pro
@mwrfits.pro
@fxaddpar.pro
@arrdelete.pro
@writefits.pro
@check_fits.pro
@sxdelpar.pro
@sxpar.pro

pro LMIRCam_remove_frame

;================================================================================

readcol, '/home/amueller/work/IDLlibs/AO/LMIRCam/datapaths.txt', tmp, format='a', /silent
print, ''
for i=0,n_elements(tmp)-1 do print, strcompress(i+1, /rem), ' ', tmp[i]

print, ''
read, 'Select Path: ', selp
path = tmp[selp-1]+'/Reduced/'

fq = ''
read, 'dx (d) / sx (s): ', fq
if (fq eq 'd') then filter = 'dx'
if (fq eq 's') then filter = 'sx'

;================================================================================

read, 'Number of bad frames: ', nbad
bad = intarr(nbad); 
read, 'Enter number of bad frames (index starts with 1!): ', bad



;PROVIDE BAD FRAMES
nbad = n_elements(bad)

bad = bad-1	;bring the numbers back to IDL counting (1st index = 0)

;================================================================================

;extract file names

file = file_search(path+'../Reduced/'+'img_'+filter+'_dc.fits', count=nfiles)
if (nfiles ne 1) then stop
pos1 = strpos(file[0], '/', /reverse_search)
file1 = strmid(file[0], pos1+1, strlen(file[0])-pos1-1)

file = file_search(path+'../Reduced/'+'vec_'+filter+'_paral.fits', count=nfiles)
if (nfiles ne 1) then stop
pos1 = strpos(file[0], '/', /reverse_search)
afile1 = strmid(file[0], pos1+1, strlen(file[0])-pos1-1)

;================================================================================

spawn, 'mv '+path+'../Reduced/'+file1+' '+path+'../Reduced/'+file1+'.orig'
spawn, 'mv '+path+'../Reduced/'+afile1+' '+path+'../Reduced/'+afile1+'.orig'

cpim1 = file1+'.orig'
acpim1 = afile1+'.orig'

;================================================================================

im1 = mrdfits(path+'../Reduced/'+cpim1, 0, hdr1, /silent)
a1 = mrdfits(path+'../Reduced/'+acpim1, 0, hdra1, /silent)
ha1 = mrdfits(path+'../Reduced/'+acpim1, 1, /silent)

sz = size(im1)

;================================================================================

good = indgen(sz[3])
for i=0,nbad-1 do good = arrdelete(good, at=bad[i]-i, length=1, /overwrite)

new1 = im1[*,*,good]

newa1 = a1[good]
newha1 = ha1[good]

; stop

; new1 = dblarr(sz[1], sz[2], sz[3]-nbad)
; new2 = dblarr(sz[1], sz[2], sz[3]-nbad)
; 
; ;assuming 1st one is OK
; if (bad[0] gt 0 and bad[1] lt sz[3]-1) then begin
; 
;   new1[*,*,0:bad[0]-1] = im1[*,*,0:bad[0]-1]
;   new1[*,*,bad[0]:*] = im1[*,*,bad[1]+1:*]
; 
;   new2[*,*,0:bad[0]-1] = im2[*,*,0:bad[0]-1]
;   new2[*,*,bad[0]:*] = im2[*,*,bad[1]+1:*]
; 
; endif else begin
; 
;   print, ''
;   print, 'Code does currently not allow the first and last frame to be bad.'
; 
; endelse

writefits, path+'../Reduced/'+file1, new1, hdr1

;================================================================================


writefits, path+'../Reduced/'+afile1, newa1, hdra1
writefits, path+'../Reduced/'+afile1, newha1, /append

;================================================================================

print, '*********************************'
print, 'Copy the new files to ../Reduced/'
print, '*********************************'

im1 = 0
new1 = 0

stop
end
