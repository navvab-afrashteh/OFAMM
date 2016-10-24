function win = CircularWindow(k)

win = zeros(2*k+1);
sze = size(win,1);

if ~(2*round(k/2)-k)
	for i = 1:sze
		for j = 1:sze
			a = pdist([i, j; k+1, k+1]);
			if (a >  min(sqrt(5*k^2/4-3*k+2),k) &&  a <= sqrt(5*k^2/4-k+1))
				win(i,j) = 1;
			end
		end
	end
	
else
	for i = 1:sze
		for j = 1:sze
			a = pdist([i, j; k+1, k+1]);
			if (a >=  min(sqrt(2)*(k-1),k) &&  a <= sqrt(5*k^2/4-k/2+1/4))
				win(i,j) = 1;
			end
		end
	end
end

