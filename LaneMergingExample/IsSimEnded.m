function IsEnded  = IsSimEnded(q, x)

safeDist = 7;

if (max(x(1), x(3)) > 20)
    IsEnded = true;
elseif (max(x(1), x(3)) > 0 && abs(x(1) - x(3)) < safeDist)
    IsEnded = true;
else
	IsEnded = false;
end

end