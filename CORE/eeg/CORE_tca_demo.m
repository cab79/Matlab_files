close all; clear; clc
addpath('C:\Data\Matlab\tensor_toolbox-master')
addpath('C:\Data\Matlab\tensor-demo-master\matlab')
pth='C:\Data\CORE\eeg\ana\prep\cleaned\part2';

select_times=[-100 499];
smooth_samples = 10;
dsample=4;
    
load('C:\Data\CORE\eeg\ana\prep\chanlocs.mat');
files = dir(fullfile(pth,'*.set'));
data=[];
for f = 1:3%length(files)
    EEG = pop_loadset('filename',files(f).name,'filepath',pth);
    data = EEG.data;
    
    dp = dsearchn(EEG.times',select_times');
    data = data(:,dp(1):dp(2),:);
    [N,T,K] = size(data)
    
    % make 2D
    data=reshape(data,N,[]);
    
    %figure 
    %hold on
    %plot(data(1:10,1:1000)','b')
    
    % smooth
    if smooth_samples
        disp('smoothing...')
        for i = 1:size(data,1)
            data(i,:) = smooth(data(i,:),smooth_samples,'moving');
        end
    end
    
    %plot(data(1:10,1:1000)','r')

    % downsample
    if dsample
        data = downsample(data',dsample)';
    end

    % make 3D
    data=reshape(data,N,[],K);
    [N,T,K] = size(data)

    %% Fit CP Tensor Decomposition

    % these commands require that you download Sandia Labs' tensor toolbox:
    % http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html

    % convert data to a tensor object
    data = tensor(data);

    % plot the ground truth
    % true_factors = ktensor(lam, A, B, C);
    % true_err = norm(full(true_factors) - data)/norm(true_factors);
    % viz_ktensor(true_factors, ... 
    %             'Plottype', {'bar', 'line', 'scatter'}, ...
    %             'Modetitles', {'neurons', 'time', 'trials'})
    % set(gcf, 'Name', 'true factors')

    % fit the cp decomposition from random initial guesses
    n_comp = 30%[4:4:64];
    err = zeros(length(n_comp),1);
    i=0;
    for n = n_comp
        R = n;
        i = i+1;
        % fit model
        est_factors = cp_als(tensor(data),R);

        % store error
        err(i) = norm(full(est_factors) - data)/norm(data);

        % visualize fit for first several fits
        if any(n == 0)%[4 8 16 32 64])
            % score aligns the cp decompositions
            %[sc, est_factors] = score(est_factors, true_factors);

            % plot the estimated factors
            viz_ktensor(est_factors, ... 
                'Plottype', {'bar', 'line', 'scatter'}, ...
                'Modetitles', {'chans', 'time', 'trials'})
            set(gcf, 'Name', ['estimated factors - fit #' num2str(n)])
            figure;
            U = est_factors.u{1};
            for u = 1:size(U,2)
                subplot(ceil(sqrt(size(U,2))),ceil(sqrt(size(U,2))),u)
                topoplot(U(:,u),chanlocs)
                title(num2str(u))
            end
        end
    end
    
    x=1;
    
end

figure(); hold on
plot(n_comp, err, 'ob')
%plot(0, true_err, 'or', 'markerfacecolor', 'r');
%xlim([-10,10])
ylim([0 1.0])
%set(gca,'xtick',[])
ylabel('model error')
legend('fits','true model')

