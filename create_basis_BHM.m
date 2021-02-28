function basis1=create_basis_BHM(N,M)
%Creates a basis of N particles with M modes

if N==1
	basis1=speye(M);
else

a=[M,ones(1,N-1)]; basis1=[];n=0;
while (a(end)<M)

	b=zeroeye(M,a(1));
	for s=2:N
		b=b+columnofones(M,a(s),a(1));
	end
	basis1=[basis1;b];
	
	a(2)=a(2)+1;
	for r=2:N-1
		if a(r)==M+1
			a(r+1)=a(r+1)+1; a(2:r)=a(r+1);
		end
	end	
	a(1)=M-a(2)+1;
	
end
basis1=[basis1;[sparse(1,M-1),N]];
end

%===========================
function mat=columnofones(M,row,len)
mat=[sparse(len,row-1),ones(len,1),sparse(len,M-row)];

%===========================
function mat1=zeroeye(M,row)
mat1=[sparse(row,M-row),speye(row)];

