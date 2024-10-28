function ans = ShowOutPuts (path)
 
 arguments 
    path="outputs/figures/"; 
 end
 msg = sprintf("input directory path, defauld is: %s\n",path); 
 x = input( msg);
 if isempty(x)
    x = path;
 end
 if isdir(path) == false
     x = path;
 end
  r = sprintf("%s/**/*.fig",x);
 listing = dir(r);
 disp(listing)
 tbl = struct2table(listing);
 tbl.date = datetime(tbl.datenum,ConvertFrom="datenum");
 tbl = removevars(tbl,"datenum");
 n = size(listing,1);
 for i=1:n
    if listing(i).isdir == false
        f = sprintf("%s%s",path,listing(i).name);
        disp(f);
        openfig(f);
        input("Enter to next ");
        close;
    end
 end

end