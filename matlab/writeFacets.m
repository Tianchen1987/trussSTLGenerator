function writeFacets(asciiMode,filename,facets)
	title='OBJECT';
	% Open the file for writing
	fid = fopen(filename,'w');

	fprintf('Writing %d bars\n',size(facets, 3));
	fprintf('...\n');
	% Write the file contents

	if asciiMode
		% Write HEADER
		fprintf(fid,'solid %s\r\n',title);
		% Write DATA
		for n=1:size(facets,3)
			fprintf(fid,[...
				'  facet normal %18.17f %18.17f %18.17f\r\n' ...
				'    outer loop\r\n' ...
				'      vertex %18.17f %18.17f %18.17f\r\n' ...
				'      vertex %18.17f %18.17f %18.17f\r\n' ...
				'      vertex %18.17f %18.17f %18.17f\r\n' ...
				'    endloop\r\n' ...
				'  endfacet\r\n'], reshape(facets(:,:,n),3,4,[]));
		end
		% Write FOOTER
		fprintf(fid,'endsolid %s\r\n',title);
	else
	% BINARY
	% Write HEADER
		binFacets=reshape(facets,12,[],size(facets,3));
		binFacets(13,:,:)=0;
		numFacets=size(binFacets,2)*size(binFacets,3);
		binFacets=reshape(binFacets,[],1,size(facets,3));
		binFacets=squeeze(binFacets);
		binFacets=reshape(binFacets,[],1);
		fprintf(fid, '%-80s', title);% Title
		fwrite(fid, numFacets, 'uint32'); % Number of facets
		% Write DATA
		for n=1:size(binFacets,1)
			if mod(n,13)~=0
				fwrite(fid, binFacets(n),'float','l');
			else
				fwrite(fid, binFacets(n),'uint16');
			end
		end
	end
	% Close the file
	fclose(fid);
	fprintf('Done writing\n');
end