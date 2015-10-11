function moyenne_dark,seq_dark
	;Computes the mean of all the dark
	;
	chemin='spectropolarimetry_m2'+'/p'+string(seq_dark)+'*.fts'
	fich=file_search(chemin,count=nf)

	for i=0,nf-1,1 do begin
			tab_dark=lecture_image(fich2(i))
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

function traitement_image,tab,tab_dark
	;Removes the dark to the raw image
	;
	result=(tab-tab_dark)
	return, result
	end

pro traitement_obs,seq_data,seq_dark
	;
	;sequence='#number of the sequences' (data,dark)
	;the .sav file contains the array with all the processed images, the sequence numbers, the dimensions of the final array and the total number of fits files
	;example:
	;traitement_obs,'002','003'
	;restore,'p002.sav'
	;

	chemin='p'+string(seq_data)+'*.fts'
	;chemin_dark='spectropolarimetry_m2'+'/p'+string(seq_dark)+'*.fts'	

	fich1=file_search(chemin,count=nf1)
	;fich2=file_search(chemin_dark,count=nf2)

	dim=size(fich1)
	dimx=dim(1)	

	result=[]
	tab_dark=moyenne_dark(seq_dark)
	
	for i=0,nf1-1,1 do begin
		tab=lecture_image(fich1(i))
		inter_result=traitement_image(tab,tab_dark)
		result=[result,inter_result]
	endfor

		
	filename='p'+string(seq_data)+'.sav'
	save,result,seq_data,seq_dark,dimx,filename=filename

	end












