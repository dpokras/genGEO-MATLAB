function [] = Assert(val1, val2, testName)

if (isnan(val1) && isnan(val2))
    disp(strcat([ testName ' ok' ]));
    return;
end
if (val1==0 && val2==0)
    disp(strcat([ testName ' ok, both 0' ]));
    return;
end

diff = abs(val1-val2)/mean([val1,val2]);
tolerance = 0.001;

if (isnan(diff))
    throw(MException('test_ResevoirDepletion:TestAssertError', strcat([ testName ' AssertError, Diff is NaN'])));
elseif (diff > tolerance)
    str = strcat([ testName ' AssertError, ' num2str(val1) ' vs ' num2str(val2) ', Diff ' num2str(diff*100) '% > ' num2str(tolerance*100) '%']);
    disp(str);
    throw(MException('test_ResevoirDepletion:TestAssertError', str));
else
    disp(strcat([ testName ' ok (' num2str(diff*100) '% diff)' ]));
end

end

