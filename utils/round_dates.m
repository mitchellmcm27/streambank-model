function [rounded] = round_dates(datevec)
%% Round dates to nearest month. Go to next month if date is later than the 15th.

year1 = year(datevec);
month1 = month(datevec);
day1 = day(datevec);

month1(day1 > 15) = month1(day1 > 15) + 1;
%month1(month1 == 13) = 1;
%year1(month1 == 13) = year1(month1 == 13) + 1; not needed for datetime

rounded = datetime(year1,month1,1);

end

