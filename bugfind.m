%Looking for bug..
clear; clc;
close all;
err_ref = importdata('out_ref.txt',' ',10);
err_ref = err_ref.data;
err = importdata('out.txt',' ',10);
err = err.data;
maxidx = 20;
subplot(2,2,1)
semilogy(1:maxidx,err_ref(1:maxidx,3),1:maxidx,err(1:maxidx,3))
%set(gca,'xdir','reverse')
subplot(2,2,2)
semilogy(1:maxidx,err_ref(1:maxidx,4),1:maxidx,err(1:maxidx,4))
%set(gca,'xdir','reverse')
subplot(2,2,3)
plot(1:maxidx,err_ref(1:maxidx,2),'b.',1:maxidx,err(1:maxidx,2),'g.')
subplot(2,2,4)
plot(1000*err_ref(:,1),err_ref(:,2),'b.',1000*err(:,1),err(:,2),'g.');
