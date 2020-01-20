function Ftot = getFtot(mol,H,X,rho,int_cord)
% [mol,H,X,info] = SCF(mol,H,X,rho,int_cord) compute force of nucleus in a dynamic 
% molecule as the gradient of Hamiltonian to find a steady state of the molecule.  
% The force will be calculated in Cartesian coordinate if zindex=0 or interior coordinate if
% zindex=1. The interior coordinate describe the relaatice pisition of 
% nucleus given by shift and rolation which do not change the structure of 
% molecule.

Fewald = getFewald(mol);
Floc = getFloc(mol,rho);
Fnl = getFnl(mol,H,X);
%compute force in Cartesian coordinate
Ftot = Fewald + Floc + Fnl;
%natom=1 case do not need  to transform coordinate
natoms=sum(mol.natoms);
if natoms<=1
    return
end
if ~exist('int_cord')
    int_cord=false;
end
if int_cord
    Ftott = Ftot';
    Ftott = Ftott(:);

    r = mol.xyzlist;
    r1 = r(1,:);
    r2 = r(2,:);
    r3 = r(3,:);
 
    %derivitave conversion of natoms=2 case
    if natoms==2
        a2 = sqrt((r(2,:)-r(1,:))*(r(2,:)-r(1,:))');
        d = (r(1,:)-r(2,:)) / a2;
        C(2,1) = [d,-d];
        fzcord=C' \ Ftott;
        Ftot=zeros(2,3);
        Ftot(2,1)=fzcord;
        return
    end

    %change xyzlist into interior coordinate, mostly can be omitted
    a2 = sqrt((r(2,:)-r(1,:))*(r(2,:)-r(1,:))');
    a3 = ((r(3,:)-r(1,:))*(r(2,:)-r(1,:))') / a2;
    b3 = sqrt((r(3,:)-r(1,:))*(r(3,:)-r(1,:))' - a3^2);
    if natoms>=4
        rr1 = repmat(r(1,:),natoms-3,1);
        a4 = ((r(4:end,:)-rr1)*(r(2,:)-r(1,:))') / a2;
        b4 = ((r(4:end,:)-rr1)*(r(3,:)-r(1,:))'  - a3*a4) / b3;
        c4 = (sqrt(diag((r(4:end,:)-rr1)*(r(4:end,:)-rr1)')' - (a4.^2)' - (b4.^2)'))';
    end

    %calculate the conversion matrix C=?z/?x, where z indicate the interior 
    %coordinate and x indicate the Cartesian coordinate
    C = zeros(3*natoms-6,3*natoms);
    C(1,1:6) = [(r1-r2)/a2 , (r2-r1)/a2];
    C(2,1:9) = [(2*r1-r2-r3)/a2 - ((r3-r1)*(r2-r1)')/a2^3*(r1-r2) , (r3-r1)/a2-a3/a2^2*(r2-r1) , (r2-r1)/a2];
    C(3,1:9) = [1/b3*((r1-r3)-a3*C(2,1:3)) , -a3*C(2,4:6)/b3 , ((r3-r1)-a3*C(2,7:9))/b3];
    if natoms>=4
        for i =4:natoms
            C(3*i-8,1:3) = ((r1-r2)+(r1-r(i,:)))/a2 - (r(i,:)-r1)*(r2-r1)'/a2^3*(r1-r2);
            C(3*i-8,4:6) = (r(i,:)-r1)/a2 - (r(i,:)-r1)*(r2-r1)'/a2^3*(r2-r1);
            C(3*i-8,3*i-2:3*i) = (r2-r1)/a2;
            C(3*i-7,1:3) = -((r3-r1)*(r(i,:)-r1)'-a3*a4(i-3))/b3^2*C(3,1:3) + 1/b3*((r1-r3)+(r1-r(i,:))-a3*C(3*i-8,1:3)-a4(i-3)*C(2,1:3));
            C(3*i-7,4:6) = -b4(i-3)/b3*C(3,4:6) - 1/b3*(a3*C(3*i-8,4:6)+a4(i-3)*C(2,4:6));
            C(3*i-7,7:9) = -b4(i-3)/b3*C(3,7:9) + 1/b3*((r(i,:)-r1)-a4(i-3)*C(2,7:9));
            C(3*i-7,3*i-2:3*i) = 1/b3*(r3-r1-a3*C(3*i-8,3*i-2:3*i));
            C(3*i-6,1:3) = 1/c4(i-3)*(r1-r(i,:)-a4(i-3)*C(3*i-8,1:3)-b4(i-3)*C(3*i-7,1:3));
            C(3*i-6,4:6) = -(a4(i-3)*C(3*i-8,4:6)+b4(i-3)*C(3*i-7,4:6))/c4(i-3);
            C(3*i-6,7:9) = -b4(i-3)*C(3*i-7,7:9)/c4(i-3);
            C(3*i-6,3*i-2:3*i) = ((r(i,:)-r1)-a4(i-3)*C(3*i-8,3*i-2:3*i)-b4(i-3)*C(3*i-7,3*i-2:3*i))/c4(i-3);
        end
    end

    %converte force Ftot into interior coordinate using C
    Ftott = C' \ Ftott;
    Ftemp = zeros(1,3*natoms);
    Ftemp(4)=Ftott(1);
    Ftemp(7)=Ftott(2);
    Ftemp(8)=Ftott(3);
    if natoms>=4
        Ftemp(10:end) = Ftott(4:end);
    end

    Ftot = reshape(Ftemp,3,natoms);
    Ftot = real(Ftot')
end

end
