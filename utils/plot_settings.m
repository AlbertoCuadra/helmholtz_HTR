function self = plot_settings()
    % Initialize struct with plot settings
    % 
    % Returns:
    %     self (struct): struct with plot settings
    
    % Description
    self.description = 'Plot settings'; 
    % Variables
    self.position = get_monitor_positions(2); % Default figure position [pixels]
    self.innerposition = [0.15 0.15 0.7 0.7]; % Set figure inner position [normalized]
    self.outerposition = [0.15 0.15 0.7 0.7]; % Set figure outer position [normalized]
    self.linestyle = '-';                     % Set line style for plots
    self.symbolstyle = 'o';                   % Set symbol style for plots
    self.linewidth = 1.8;                     % Set line width for plots
    self.fontsize = 20;                       % Set fontsize
    self.colorpalette = 'Seaborn';            % Set Color palette (see brewermap function for more options)
    self.colorpaletteLenght = 11;             % Set Maximum number of colors to use in the color palette
    self.box = 'off';                         % Display axes outline
    self.grid = 'off';                        % Display or hide axes grid lines
    self.hold = 'on';                         % Retain current plot when adding new plots
    self.axis_x = 'tight';                    % Set x-axis limits
    self.axis_y = 'auto';                     % Set y-axis limits
    self.xscale = 'linear';                   % Set x-axis scale (linear or logarithmic)
    self.yscale = 'linear';                   % Set y-axis scale (linear or logarithmic)
    self.xdir = 'normal';                     % Set x-axis direction (normal or reverse)
    self.ydir = 'normal';                     % Set y-axis direction (normal or reverse)
    self.title = [];                          % Set title
    self.label_type = 'medium';               % Set label with variable (short), name (medium), or name and variable (long)
    self.labelx = [];                         % Set x label
    self.labely = [];                         % Set y label
    self.legend_name = [];                    % Set legend labels
    self.legend_location = 'northeastoutside';% Set legend location
    self.colorline = [44, 137, 160]/255;      % Default colorline
    self.colorlines = [135, 205, 222;...      % Default colorlines
                       95, 188, 211;...       % Default colorlines
                       44, 137, 160;...       % Default colorlines
                       22,  68,  80]/255;     % Default colorlines
    self.blue = [0.3725, 0.7373, 0.8275];     % Default colorline
    self.gray = [0.50, 0.50, 0.50];           % Default colorline
    self.red = [0.64,0.08,0.18];              % Default colorline
    self.orange = [212, 85, 0]/255;           % Default colorline
    self.brown = [200, 190, 183]/255;         % Default colorline
    self.brown2 = [72, 55, 55]/255;           % Default colorline
end