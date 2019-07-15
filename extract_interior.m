function [t_interior,f_interior] = extract_interior(t_data,f_data,t_min,t_max)
%EXTRACTS all points in range tmin<= t_data <= tmax, and corresponding
%points of f (of same length)
    Ndata=length(t_data);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if t_data(i)>=t_min
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if t_data(i)<=t_max
            Nmax=i;
        end
    end
    t_interior=t_data(1,Nmin:Nmax);
    f_interior=f_data(1,Nmin:Nmax);
end

