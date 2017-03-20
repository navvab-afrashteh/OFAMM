function [CorrContour] = contourConfirmation(C, sVerified, rs, cs, Nnested)

CorrContour = sVerified;
for n = size(sVerified,2):-1:1
    k = sVerified(1,n);
    N = sVerified(2,n);
    c = C(:,k+1:k+N);
    x = c(1,:);
    y = c(2,:);
    [inn,~] = inpolygon(cs,rs, x,y);
    if ~inn
        CorrContour(:,n) = [];
    end
end

% Refining Sources
if size(CorrContour,2) >= Nnested
	sVerified = CorrContour(:,1);
	
	for t = 2:size(CorrContour,2)
		k1 = CorrContour(1,t);
		n1 = CorrContour(2,t);
		c1 = C(:, k1+1:k1+n1);
		
		q=0;
		len = size(sVerified,2);
		for p = 1:len
			k2 = sVerified(1,p);
			n2 = sVerified(2,p);
			c2 = C(:, k2+1:k2+n2);
			
			[inside, Cin] = isInsideContour(c1 , c2);
			if inside == 1
				if size(Cin,2) == size(c1,2)
					d = Cin - c1;
					d = d(:);
					if ~sum(d)
						sVerified(:,p) = [];
						sVerified = [sVerified, CorrContour(:,t)];
					end
				end
				break
				
			elseif inside == 0
				q = q+1;
			end
		end
		
		if q == len
			sVerified = [sVerified, CorrContour(:,t)];
		end
	end
else
	sVerified = [];
end
CorrContour = sVerified;