function updateplot(hfig,mol,H,isoval,alphaval,fatoms,fbonds, ...
    fdensity,fcenterdensity)

if fatoms
    vatom = 'on';
else
    vatom = 'off';
end

if fbonds
    vbond = 'on';
else
    vbond = 'off';
end

for it = 1:length(hfig.Children)
    if strcmpi(hfig.Children(it).UserData,'bond')
        hfig.Children(it).Visible = vbond;
        continue;
    end
    if strcmpi(hfig.Children(it).UserData,'atom')
        hfig.Children(it).Visible = vatom;
        continue;
    end
    if strcmpi(hfig.Children(it).UserData,'density')
        
        if fdensity
            delete(hfig.Children(it));
            plotdensity(hfig,mol,H,isoval,alphaval,fcenterdensity);
        else
            hfig.Children(it).Visible = 'off';
        end
        continue;
    end
end

end