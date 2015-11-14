function moyenne_dark,seq_dark
	;Computes the mean of all the dark
	;
	chemin='spectropolarimetry_m2'+'/p'+string(seq_dark)+'*.fts'
	fich=file_search(chemin,count=nf)

	for i=0,nf-1,1 do begin
			tab_dark=lecture_image(fich(i))
			dark=[dark,tab_dark]
	endfor

	dim=size(dark)
	dimx=dim(1)
	dimy=dim(2)
	mean_dark=make_array(688,520)
	for j=0,dimy-1,1 do begin
		for i=0,687,1 do begin
			x=[]
			for k=0,nf-1,1 do begin
				x=[x,dark(i+k*688,j)]
			endfor
			mean_dark(i,j)=mean(x)
		endfor
	endfor
	return,mean_dark
	end

function moyenne_flat,seq_flat,seq_flat_dark
	;Computes the mean of all the dark
	;
	chemin1='spectropolarimetry_m2'+'/p'+string(seq_flat)+'*.fts'
	fich1=file_search(chemin1,count=nf1)
	chemin2='spectropolarimetry_m2'+'/p'+string(seq_flat_dark)+'*.fts'
	fich2=file_search(chemin2,count=nf2)

	for i=0,nf1-1,1 do begin
			tab_flat=lecture_image(fich1(i))
			tab_flat_dark=lecture_image(fich2(0))
			flat=[flat,tab_flat-tab_flat_dark]
	endfor

	dim=size(flat)
	dimx=dim(1)
	dimy=dim(2)
	mean_flat=make_array(688,520)
	for j=0,dimy-1,1 do begin
		for i=0,687,1 do begin
			x=[]
			for k=0,nf-1,1 do begin
				x=[x,flat(i+k*688,j)]
			endfor
			mean_flat(i,j)=mean(x)
		endfor
	endfor
	return,mean_flat
	end


function traitement_image,tab,tab_dark
	;Removes the dark to the raw image (specific to way changing)
	;
	result=(tab-tab_dark)
	return, result
	end

function traitement_image_bis,tab,tab_dark,tab_flat,tab_flat_dark
	;Usual treatment on a raw image
	;
	result=(tab-tab_dark)/(tab_flat-tab_flat_dark)
	return, result
	end

function lecture_image,filename,result
	tab=readfits(filename,header)
	tab=uint(tab)
	tab=float(tab)
	return, tab
	end

pro traitement_obs,seq_data,seq_dark,seq_flat,seq_flat_dark
	;
	;sequence='#number of the sequences' (data,dark,flat,dark_flat)
	;the .sav file contains the array with all the processed images, the sequence numbers, the dimensions of the final array and the total number of fits files
	;example:
	;traitement_obs,'002','003','004','005'
	;restore,'p002.sav'
	;

	chemin='p'+string(seq_data)+'*.fts'
	chemin_flat='p'+string(seq_flat)+'*.fts'	

	fich1=file_search(chemin,count=nf1)
	fich2=file_search(chemin_flat,count=nf2)

	dim=size(fich1)
	dimx=dim(1)	

	result=[]
	tab_dark=moyenne_dark(seq_dark)
	tab_flat=moyenne_flat(seq_flat,seq_flat_dark)
	
	if string(seq_data) EQ '001' then begin
		for i=0,nf1-1,1 do begin
			tab=lecture_image(fich1(i))
			inter_result=traitement_image(tab,tab_dark)
			result=[result,inter_result]
		endfor
	endif else begin
		for i=0,nf1-1,1 do begin
			tab=lecture_image(fich1(i))
			flat=lecture_image(fich2(0))
			inter_result=traitement_image_bis(tab,tab_dark,flat,flat_dark)
			result=[result,inter_result]
		endfor
	endelse

	intermediate=result(0:688,*)
	window,0,xs=688,ys=520
	tvscl,intermediate
	print,'Cliquez sur la première ligne de référence'
	cursor,ref_x0,ref_y0
	;real_y0=where(min(tab(ref_x0,ref_y0-2:ref_y0+2)))
	print,'Cliquez sur la deuxième ligne de référence'
	cursor,ref_x1,ref_y1
	;real_y1=where(min(tab(ref_x1,ref_y1-2:ref_y1+2)))
	height=ref_y0-ref_y1
	print,'Hauteur de chaque ligne : ',height
	print,'Cliquez sur la ligne (la plus brillante du couple) pour laquelle le profil sera tracé ainsi que celui de la ligne associée'
	cursor,line1_x,line1_y
	line2_y=line1_y+height
	;print,'Exemple de tracé de profils'
	prof1=tab(*,line1_y)
	prof2=tab(*,line2_y)
	;prof3=(prof1+1.94*prof2)/2.
	;window,1,xs=688,ys=300
	cgplot,prof1,ytitle='ADU',xtitle='Pixels',title='Exemple de profils'
	oplot,1.94*prof2,color=cgcolor('red')
	;oplot,prof3,color=cgcolor('blue')
	dim=size(prof1)
	print,dim(1),dim(2)
	
	profil1=[]
	profil2=[]
	prof1=0.
	prof2=0.
	for i=0,(nf1-2)*687.,687 do begin
		;print,i,i+1
		tab=result(i:(i+687),*)
		prof1=tab(*,line1_y)
		prof2=tab(*,line2_y)
		profil1=[profil1,prof1]
		profil2=[profil2,prof2]
	endfor
	;tab=result((nf1-1)*687.:*,*)
	;profil1=[profil1,tab(*,line1_y)]
	;profil2=[profil2,tab(*,line2_y)]



	filename='p'+string(seq_data)+'.sav'
	save,result,seq_data,seq_dark,dimx,height,profil1,profil2,filename=filename

	end

