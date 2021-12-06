;---read topo
ftopo=addfile("./DATA/GTOPO30_10MIN.CDF","r")
hh0=ftopo->HT(0,:,:)
hh=short2flt(hh0)
copy_VarCoords(hh0, hh)

;---read files
mis	= -999.
fin1	= "./temp_info.txt"
fin2    = "./temp.txt"

readin	 = asciiread(fin1,(/3/),"string")

year    = str_get_cols(readin(0),1,4)
imonth	= stringtofloat(str_get_cols(readin(0),5,6))
iday	= stringtofloat(str_get_cols(readin(0),7,8))
ihour	= stringtofloat(str_get_cols(readin(0),9,10))
iminute	= stringtofloat(str_get_cols(readin(0),11,12))
if(imonth.lt.10) then
	 month	= "0"+flt2string(imonth)
else
	 month	= flt2string(imonth)
end if
if(iday.lt.10) then
	 day	= "0"+flt2string(iday)
else
	 day	= flt2string(iday)
end if
if(ihour.lt.10) then
	 hour	= "0"+flt2string(ihour)
else
	 hour	= flt2string(ihour)
end if
if(iminute.lt.10) then
	 minute	= "0"+flt2string(iminute)
else
	 minute	= flt2string(iminute)
end if

date    = year+month+day+hour
lonW	= 	stringtofloat(str_get_cols(readin(1),0,14))
lonE	= stringtofloat(str_get_cols(readin(1),15,29))
latS	= 	stringtofloat(str_get_cols(readin(1),30,44))
latN	= stringtofloat(str_get_cols(readin(1),45,59))
NX	= 	stringtoint(str_get_cols(readin(2),0,14))
NY	= stringtoint(str_get_cols(readin(2),15,29))
NT	= stringtoint(str_get_cols(readin(2),30,35))

grid   = asciiread(fin2,(/NT,NY,NX/),"float")
lat	= fspan(latS,latN,NY)
lon	= fspan(lonW,lonE,NX)

lat!0  = "lat"
lat@long_name  = "lat"
lat@units      = "degrees_north"
lat@lat        = lat
lon!0  = "lon"
lon@long_name  = "lon"
lon@units      = "degrees_east"
lon@lon        = lon
time	= ispan(0,(NT-1)*6,6)
time@calendar	= "standard"
time@units	= "minutes since "+year+"-"+month+"-"+day+" "+hour+":"+minute+":00"
time@standard_name	= "time"
utc	= cd_calendar(time,0)

utc_date=new(NT,"string")
utc_date(:)=sprinti("%0.4i",tointeger(utc(:,0)))+sprinti("%0.2i",tointeger(utc(:,1)))+ \
						sprinti("%0.2i",tointeger(utc(:,2)))+sprinti("%0.2i",tointeger(utc(:,3)))+ \
						sprinti("%0.2i",tointeger(utc(:,4)+utc(:,5)/60))


grid@_FillValue	= mis
grid!0	= "time"
grid!1	= "lat"
grid!2	= "lon"
grid&time	= time
grid&lat	= lat
grid&lon	= lon

;---write into nc file
root	= "./DATA/"
fout    = root+"MCS_grid_FY2E_Xinjiang_"+date+".nc"
print("fout:"+fout)
system("rm -f "+fout)
ncdf	= addfile(fout, "c")

;---define dimension
fAtt= True
fAtt@title= "MCS FY2E grid data"
fAtt@source_file= "FY2E_TBB_AWX"
fAtt@Conventions= "None"
fAtt@creation_date= systemfunc("date")
fileattdef(ncdf,fAtt)

filedimdef(ncdf,"time",-1,True)

ncdf->grid	= grid


;---plot section
grid=where(grid.eq.0, grid@_FillValue, grid)

mpres=True
mpres@gsnFrame=False
mpres@gsnDraw=False
mpres@gsnAddCyclic=False
mpres@cnFillOn=True
mpres@cnLinesOn=False
mpres@cnLineLabelsOn=False
mpres@cnInfoLabelOn=False
mpres@mpMinLatF=latS
mpres@mpMaxLatF=latN
mpres@mpMinLonF=lonW
mpres@mpMaxLonF=lonE
mpres@lbLabelBarOn=False
mpres@cnLevelSelectionMode="ExplicitLevels"
mpres@cnLevels=ispan(500,5000,500)
mpres@cnFillPalette="gsltod"
mpres@gsnRightString=year+"."+month+"."+day+" "+hour+":00"


res=True
res@gsnFrame=False
res@gsnDraw=False
res@gsnAddCyclic=False
res@cnFillOn=True
res@cnLinesOn=False
res@cnLineLabelsOn=False
res@cnInfoLabelOn=False
res@cnConstFEnableFill=True
res@cnConstFLabelOn=False
res@cnLevelSelectionMode="ExplicitLevels"
res@cnFillPalette="MPL_Spectral"
res@cnFillOpacityF=0.5
res@lbLabelBarOn=False
if(max(grid).gt.1) then
  	res@cnLevels=ispan(0,toint(max(grid))-1,1)+0.5
else
	res@cnLevels=(/0.5/)
end if
wks=gsn_open_wks("eps", "MCS_"+date)
ploth=gsn_csm_contour_map(wks, hh({latS-0.1:latN+0.1},{lonW-0.1:lonE+0.1}), mpres)
plot=gsn_csm_contour(wks, grid(0,:,:), res)
overlay(ploth, plot)
draw(ploth)
frame(wks)