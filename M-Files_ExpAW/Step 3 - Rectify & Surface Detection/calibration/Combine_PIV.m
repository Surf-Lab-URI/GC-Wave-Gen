function [ComboPIV] = Combine_PIV(PIV_A,PIV_W, Ix, PF_Surface)

Size1 = size(PIV_W,1); %max(size(PIV_A,1),size(PIV_W,1));
Size2 = size(PIV_W,2); %max(size(PIV_A,2),size(PIV_W,2));

ComboPIV = NaN(Size1,Size2);
Lx = size(PIV_A,2);
PIV_A = interp2((1:Lx),(1:size(PIV_A,1))',PIV_A,Ix,(1:size(PIV_A,1))');
%PIV_A(:,[Ix<0,Ix>Lx]) = NaN;

for i = 1:size(ComboPIV,2)
    if i <= size(PIV_A,2)
        ComboPIV(1:floor(PF_Surface(i)),i) = PIV_A(1:floor(PF_Surface(i)),i);
    end
    ComboPIV(floor(PF_Surface(i))+1:end,i) = PIV_W(floor(PF_Surface(i))+1:end,i);
end

%For PIVLAB ad lines below
%Mask(Mask==1)=0;
%Mask(isnan(Mask))=1;

end