clear all;

angle = 360*rand(9,3);

aeuler = zeros(81,3);

for i=1:9
	for j=1:9
		id=9*(i-1) + j;

		x = ceil(i/3);
		y = ceil(j/3);
		grain = 3*(x-1)+y;
		aeuler(id,:) = angle(grain,:);
	end
end

save aeuler aeuler -ascii
