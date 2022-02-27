function [facets,listNode,listElem]=genFacets(x1,x2,div,r)
	if norm(x1)>norm(x2)
		temp=x1;
		x1=x2;
		x2=temp;
	end
	d=x2-x1;
	projLen=sqrt(d(1)^2+d(3)^2);
	phi=atan(d(2)/projLen);
	theta=acos(d(3)/projLen);
	if d(1)<0
		theta=2*pi-theta;
	end
	
	if projLen~=0
		M1=makehgtform('yrotate',theta);
	else
		M1=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
	end
	x_p_axis=M1*[1 0 0 0]';
	M2=makehgtform('axisrotate',-x_p_axis(1:3),phi);
	
	M3=makehgtform('translate',x1);
	M4=makehgtform('translate',x2);
	
	M_1=M3*M2*M1;
	M_2=M4*M2*M1;
	
	angle=0:pi/div:2*pi-pi/div;
	areaPoly=r^2*pi;
	sPoly=sqrt(4*areaPoly/(div*2)/cot(pi/(div*2)));
	rPoly=sPoly*csc(pi/(div*2))/2;
	
	x=rPoly*cos(angle);
	y=rPoly*sin(angle);
	
	x_coord=[x' y' zeros(length(x),1) ones(length(x),1)];
	
	for n=1:length(x)
		x_trans_1(n,:)=M_1*x_coord(n,:)';
		x_trans_2(n,:)=M_2*x_coord(n,:)';
	end
    
	listNode=[x_trans_1;x_trans_2];
	num.Elem=div*2;

	n_start=num.Elem*2+1;
	n_end=num.Elem*2+2;
	
	listNode(n_start,:)=[x1,1];
	listNode(n_end,:)=[x2,1];

	eC=1;
	listElem=zeros(num.Elem,3);
	for n=1:num.Elem
		if n~=num.Elem
			n1=n;
			n2=n+1;
			n3=n+num.Elem;
			n4=n+num.Elem+1;
			listElem(eC,:)=[n1 n2 n3];
			eC=eC+1;
			listElem(eC,:)=[n2 n4 n3];
			eC=eC+1;
			listElem(eC,:)=[n2 n1 n_start];
			eC=eC+1;
			listElem(eC,:)=[n3 n4 n_end];
			eC=eC+1;
		else
			n1=n;
			n2=n-num.Elem+1;
			n3=n+num.Elem;
			n4=n+1;
			listElem(eC,:)=[n1 n2 n3];
			eC=eC+1;
			listElem(eC,:)=[n2 n4 n3];
			eC=eC+1;
			listElem(eC,:)=[n2 n1 n_start];
			eC=eC+1;
			listElem(eC,:)=[n3 n4 n_end];
			eC=eC+1;
		end
	end

	facets=zeros(3,4,num.Elem*4);
	for n=1:num.Elem*4
		
		v1=listNode(listElem(n,2),1:3)-listNode(listElem(n,1),1:3);
		v2=listNode(listElem(n,3),1:3)-listNode(listElem(n,1),1:3);
		vNorm=cross(v1,v2);
		vNorm=vNorm/norm(vNorm);
		temp=[
			vNorm;...
			listNode(listElem(n,1),1:3);...
			listNode(listElem(n,2),1:3);...
			listNode(listElem(n,3),1:3);
		];
		facets(:,:,n)=temp';
	end
	facets=reshape(facets,3,[],1);
end
