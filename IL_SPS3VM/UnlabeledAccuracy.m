accuracyBL_t = [0.3261,0.4980,0.1141,0.3583,0.6012,0.5025,0.4787,0.2564];
accuracyBC_t = [0.3311,0.4961,0.0980,0.4327,0.5902,0.9089,0.7541,0.8450];
accuracyIL_BC_t = [0.3311,0.4963,0.0980,0.5016,0.5850,0.9109,0.6783,0.8462];
accuracyBL = [];
accuracyBC = [];
accuracyIL_BC= [];
for i=1:10
    accuracyBL(i,:) = accuracyBL_t;
    accuracyBC(i,:) = accuracyBC_t;
    accuracyIL_BC(i,:) = accuracyIL_BC_t;
end


x = [1,2,3,4,5,6,7,8];
hFig = figure(1);
set(hFig, 'Position',  1.15*[232   206   560   380])    
    
 	%boxplotCsub( dataset(:,:,3),1,'+',1,1,'b',0,2,true,[3 3])
   %axis;axis([min(min(min(dataset))) max(max(max(dataset))) ans(3) ans(4)])


ylabel('Unlabeled Accuracy','fontsize',14);
p = plot(x,accuracyBL_t,'+k',x,accuracyBC_t,'+b',x,accuracyIL_BC_t,'*r')
set(p,'Linewidth',2)

h = legend('BL-S^3VM','BCS^3VM', 'IL-BCS^3VM');


set(h,'Location','northwest')
set(h,'String',{'BL-S^3VM','BCS^3VM', 'IL-BCS^3VM'}) 
    
xlabel('','fontsize',13);
set(gca,'xticklabel',{'CodRNA','Madelon','IJCNN1','Text','Usps','W6a','A9a','Mushrooms'},'fontsize',11);
axis([0.5 8.5 0 1.00]);
filename=['F:\nkp\ILS3VM_blance_xm\','UnlabeledAccuracy.jpg'];
saveas(hFig,filename);