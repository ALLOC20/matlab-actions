clc;
clear;
close all;

for ii=1:503
    tic
    fprintf('index = %d,', ii);
    if ii<10
        index="000"+ii;
    else
        if (ii>9)&&(ii<100)
            index="00"+ii;
        else
            if ii>99
                index="0"+ii;
            end
        end
    end
    im_name1 = "image/out_of_focus" + index + ".jpg";
    Iuint1 = imread(im_name1);
    im_name2 = "image/out_of_focus" + index + ".png";
    Iunit2= imread(im_name2);
    Idouble1 = im2double(Iuint1);
    Idouble2 = im2double(Iunit2);
    [psf, sparse_psf, reliable_edge_map] = blur_estimate_our(Idouble1);
    psf_1=psf/max(max(psf));
    [m,n]=size(psf_1);
    
    
    u=0.5;

    J=zeros(m,n);
    for i=1:m
        for j=1:n
            if psf_1(i,j)<u
                J(i,j)=0;
            else
                J(i,j)=1;
            end
        end
    end
    K=zeros(m,n);
    for i=1:m
        for j=1:n
            if (Idouble2(i,j)==J(i,j))&&(J(i,j)==0)
                K(i,j)=0;
            else
                K(i,j)=1;
            end
        end
    end
%     figure(1);
%     subplot(2,2,1), imshow(Idouble1); title('input');
%     subplot(2,2,2), imshow(J); title('mypsf');
%     subplot(2,2,3), imshow(Idouble2); title('given');
%     subplot(2,2,4), imshow(K); title('and');
    R=sum(sum(sum(J==0)));
    Rg=sum(sum(sum(Idouble2==0)));
    T=sum(sum(sum(K==0)));
    p(ii)=T/R;
    r(ii)=T/Rg;
    figure(2);
    plot(r(ii),p(ii),'ro');
    hold on;
    toc
end

figure(3);
axis([0 1 0.5 1]);
hold on
plot(r,p,'-o');
xlabel('recall');
ylabel('precision');

save data.mat p r;
print('-f2','-dpng','pic2.png');
print('-f3','-dpng','pic3.png');



