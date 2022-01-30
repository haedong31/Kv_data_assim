function seq = spgetseq(n,d)
% SPGETSEQ internal function that returns all a set of indices in the
% natural numbers and 0 such that they are less than n
%Code follows from Andreas Klimke's  "Sparse Grid Interpolation Toolbox"
%v1.1
nlevels = nchoosek(double(n+d-1),double(d-1));
seq = zeros(nlevels,d);
seq(1,1) = n;
max = n;
for k = 2:nlevels
	if seq(k-1,1) > 0
		seq(k,1) = seq(k-1,1) - 1;
		for l = 2:d
			if seq(k-1,l) < max
				seq(k,l) = seq(k-1,l) + 1;
				for m = l+1:d
					seq(k,m) = seq(k-1,m);
				end
				break;
			end
		end
	else
		sum = 0;
		for l = 2:d
			if seq(k-1,l) < max
				seq(k,l) = seq(k-1,l) + 1;
				sum = sum + seq(k,l);
				for m = l+1:d
					seq(k,m) = seq(k-1,m);
					sum = sum + seq(k,m);
				end
				break;
			else
				temp = 0;
				for m = l+2:d
					temp = temp + seq(k-1,m);
				end
				max = n-temp;
				seq(k,l) = 0;
			end
		end
		seq(k,1) = n - sum; max = n - sum;
	end
end
