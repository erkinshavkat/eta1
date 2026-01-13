dat20Hz = load('results/highlowfreq/20Hz.mat');
dat40Hz = load('results/highlowfreq/40Hz.mat');
dat80Hz = load('results/highlowfreq/80Hz.mat');
dat160Hz = load('results/highlowfreq/160Hz.mat');
dat240Hz = load('results/highlowfreq/240Hz.mat');
dat480Hz = load('results/highlowfreq/480Hz.mat');
dat960Hz = load('results/highlowfreq/960Hz.mat');

datas= {dat20Hz, dat40Hz, dat80Hz, dat160Hz, dat240Hz, dat480Hz, dat960Hz};
freqs = [20,40,80,160,240,480,960];

v = VideoWriter('vis/H A14 var omega_0.mp4','MPEG-4');
v.FrameRate = 6;
open(v);


nimpacts=50;
fig = figure('Position', [0, 0, 1500, 900]);
for n=1:nimpacts
    clf;
    hold on;
    for i=1:length(datas)
        data = datas{i};
        freq = freqs(i);
        eta_data = data.eta_data;
        H_data = data.H_data;
        Nx = data.Nx; Lx = data.Lx; Nk = data.Nk; Lk = data.Lk; hx = Lx/Nx;
        x = hx*(0:Nx-1)-Lx/2;
        K_vec = linspace(0,Lk,Nk)';
        dk= K_vec(2)-K_vec(1);
        K_vec=K_vec+dk;

        fft_domain = (0:Nx-1)/Lx;

        H_full_ax = plot(K_vec/(2*pi), H_data(:,n),  'LineWidth', 2);  
        xlim([0,5])
        ylim([-0.2 0.2])
        xlabel('x/\lambda_F'); ylabel('\eta');


        % H_full_ax = plot(x, eta_data(:,n),  'LineWidth', 2);  
        % xlim([-5,5])
        % ylim([-0.01 0.01])
        % xlabel('x/\lambda_F'); ylabel('\eta');

        % fftspectrum= abs(fft(eta_data(:,n)));
        % H_full_ax = plot(fft_domain, fftspectrum,  'LineWidth', 2);  
        % xlim([0,3])
        % ylim([0,0.02])
        % xlabel('k/k_F'); ylabel('FFT(|\eta|)');


        legend_strings{i} = sprintf('%d Hz', freq);
    end
    legend(legend_strings);
    title(sprintf('Impact number: %d', n));
    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);