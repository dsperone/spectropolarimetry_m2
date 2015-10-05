pro traitement_obs,seq_data,seq_dark;,seq_flat,seq_dark_flat
	;sequence='#number of the sequences' (data,dark,flat,falt_dark)
	;example:
	;traitement_obs,'002'
	chemin='spectropolarimetry_m2'+'/p'+string(seq_data)+'*.fts'
	chemin_dark='spectropolarimetry_m2'+'/p'+string(seq_dark)+'*.fts'	
	;chemin_flat='spectropolarimetry_m2'+'/p'+string(seq_flat)+'*.fts'
	;chemin_flat_dark='spectropolarimetry_m2'+'/p'+string(seq_dark_flat)+'*.fts'

	fich1=file_search(chemin,count=nf1)
	fich2=file_search(chemin_dark,count=nf2)
	;fich3=file_search(chemin_flat,count=nf3)
	;fich4=file_search(chemin_flat_dark,count=nf4)

	;if (sequence eq '01') OR (sequence eq '02') then begin
	;	chemin_dark='spectropolarimetry_m2'+'/p003*x1.fts'
	;	chemin_flat='spectropolarimetry_m2'+'/p004*y1.fts'
	;	chemin_flat_dark='spectropolarimetry_m2'+'/p005*x1.fts'
	;endif else if (sequence eq '11') then begin
	;	chemin_dark='spectropolarimetry_m2'+'/p012*x1.fts'
	;	chemin_flat='spectropolarimetry_m2'+'/p013*y1.fts'
	;	chemin_flat_dark='spectropolarimetry_m2'+'/p012*x1.fts'
	;endif else begin
	;	chemin_dark='spectropolarimetry_m2'+'/p022*x1.fts'
	;	chemin_flat='spectropolarimetry_m2'+'/p013*y1.fts'
	;	chemin_flat_dark='spectropolarimetry_m2'+'/p022*x1.fts'
	;endelse

	result=[]

	tab_dark=lecture_image(fich2(0))
	;tab_flat=lecture_image(fich3(0))
	;tab_flat_dark=lecture_image(fich4(0))
	dim=size(fich1)
	dimx=dim(1)	
	
	for i=0,dimx-1,1 do begin
		tab=lecture_image(fich1(i))
		inter_result=traitement_image(tab,tab_dark)
		result=[result,inter_result]
	endfor

		
	filename='p'+string(seq_data)+'.sav'
	save,result,seq_data,seq_dark,filename=filename

	end












