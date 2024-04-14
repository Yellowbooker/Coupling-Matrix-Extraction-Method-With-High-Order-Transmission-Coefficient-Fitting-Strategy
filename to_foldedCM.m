function foldedCM=to_foldedCM(N,M)

% It is a function to tranfrom the transversal matrix to the folded matrix.
% This code including the normalized operation to make sure all of the 
% main-coupling are positive.
% By YB 
% Date: 2024-04-14

    N = N+2;
    for i=1:N-3 
        if mod(i,2)==1
            c=-1;
            for j=1:N-3+1-i
                l=N-(i+1)/2-j+1;
                nn=l-1;
                k=(i+1)/2;
                m=k;
                theta=atan(c*M(k,l)/M(m,nn));
                R=eye(N);
                R(nn,nn)=cos(theta);
                R(l,l)=R(nn,nn);
                R(nn,l)=-sin(theta);
                R(l,nn)=-R(nn,l);
                M=R*M*R.';
            end
        else
            c=1;
            for j=1:N-3+1-i
                k=3+i/2+j-2;
                m=k+1;
                l=N-i/2+1;
                nn=l;
                theta=atan(c*M(k,l)/M(m,nn));
                R=eye(N);
                R(k,k)=cos(theta);
                R(m,m)=R(k,k);
                R(k,m)=-sin(theta);
                R(m,k)=-R(k,m);
                M=R*M*R.';
            end
        end
    end
    % normalized operation
    for i = 1:N-1
        if real(M(i,i+1)) < 0
            R = eye(N,N);
            R(i+1,i+1) = -1;
            M = R*M*R.';
        end
    end
    for i = 1:N
        for j = 1:N
            if abs(M(i,j)) < 1e-4
                M(i,j) = 0;
            end
        end
    end

    foldedCM = M;
    
end
