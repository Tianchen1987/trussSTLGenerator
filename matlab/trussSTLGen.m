function trussSTLGen(fPath,fCoord,fCon,fOut,div, asciiMode)

	%	input:
	%	fPath: path of the input file
	%	fCoord: file name of the coordinate file
	%	fCon: file name of the connectivity file
	%	fOutL file name of the STL
	%	iDiv: number of polygon to approximate the circular cross section, 4 would mean a square
	%	boolAscii: True - ascii, false - binary (smaller file size but not human readable)
	
	sN=csvread(strcat(fPath,fCoord)); % id, x, y, z
	sE=csvread(strcat(fPath,fCon)); % id, n1, n2

	numBars=size(sE,1);

	facetsObj=zeros(3,div*32,numBars);

	for n=1:numBars
		x1=[sN(sE(n,1),1) sN(sE(n,1),2) sN(sE(n,1),3)];
		x2=[sN(sE(n,2),1) sN(sE(n,2),2) sN(sE(n,2),3)];
		[facetsObj(:,:,n),~,~]=genFacets(x1,x2,div,sE(n,3));
	end

	writeFacets(asciiMode,strcat(fPath,fOut),facetsObj);
end