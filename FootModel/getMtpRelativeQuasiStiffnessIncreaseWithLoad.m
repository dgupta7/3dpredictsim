
function [k_mtj, mtprqsiwl] = getMtpRelativeQuasiStiffnessIncreaseWithLoad(R)

BW = R.S.mass*9.81;
n_tib = length(R.Fs_tib);

for i=1:n_tib
        js = find(R.failed(:,i)==0);
        T_mtp = R.M_mtp(js,i);
        
        pl = polyfit(R.Qs_mtp(js)'*180/pi,T_mtp,1);
        pln = polyfit(R.Qs(js,i,R.jointfi.mtp.r)*180/pi,T_mtp,1);

        kg(i) = -pl(1);
        kf(i) = -pln(1);

        if isfield(R,'T_mtp')
            ple = polyfit(R.Qs_mtp(js)'*180/pi,R.T_mtp(js,i),1);
            plne = polyfit(R.Qs(js,i,R.jointfi.mtp.r)*180/pi,R.T_mtp(js,i),1);

            kge(i) = ple(1);
            kfe(i) = plne(1);
        end
end
    

pln = polyfit(R.Fs_tib/BW,kfe./kfe(1),1);
disp(num2str(pln(1),5))
        

 mtprqsiwl = pln(1);
 k_mtj = R.S.kMT_li;















