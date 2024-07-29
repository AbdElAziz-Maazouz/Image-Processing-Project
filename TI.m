function varargout = TI(varargin)
% TI MATLAB code for TI.fig
%      TI, by itself, creates a new TI or raises the existing
%      singleton*.
%
%      H = TI returns the handle to a new TI or the handle to
%      the existing singleton*.
%
%      TI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TI.M with the given input arguments.
%
%      TI('Property','Value',...) creates a new TI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TI

% Last Modified by GUIDE v2.5 25-Feb-2024 12:48:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TI_OpeningFcn, ...
                   'gui_OutputFcn',  @TI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before TI is made visible.
function TI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TI (see VARARGIN)

% Choose default command line output for TI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Menu_Fichier_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to Ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [file, path] = uigetfile('*.*');
    handles.ima = imread(fullfile(path, file)); % Use fullfile to construct full path
    handles.courant_data = handles.ima;

    % Display the original image in imgo axes
    axes(handles.imgO);
    imshow(handles.courant_data);

    handles.output = hObject;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Enregister_Callback(hObject, eventdata, handles)
% hObject    handle to Enregister (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    image = handles.ima_traite;
    [file,path] = uiputfile('*.png','Enregistrer Votre Image ... ');
    imwrite(image,sprintf('%s',path,file),'png');

% --------------------------------------------------------------------
function Quitter_Callback(hObject, eventdata, handles)
% hObject    handle to Quitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    delete(handles.figure1)


% --------------------------------------------------------------------
function Menu_About_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Créer le message à afficher
    message = {'Ce projet réalisé par : Maazouz Abd El Aziz', ...
               'Encadré par : M. Hamid Tahiri', ...
               'Master : MLAIM 2023/2024'};

    % Afficher la fenêtre modale avec le message
    msgbox(message, 'À propos', 'modal');


% --------------------------------------------------------------------
function FiltrageSpatial_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrageSpatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltrageFrequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrageFrequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Morphologie_Callback(hObject, eventdata, handles)
% hObject    handle to Morphologie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PointsDinteret_Callback(hObject, eventdata, handles)
% hObject    handle to PointsDinteret (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Transformation_Callback(hObject, eventdata, handles)
% hObject    handle to Transformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Binarisation_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Appliquer la binarisation en utilisant un seuil (ajustez selon vos besoins)
    seuil = 128;
    image_binarisee = image > seuil;

    % Convertir l'image binaire en uint8 pour l'affichage
    image_binarisee = uint8(image_binarisee * 255);

    % Afficher l'image binarisée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_binarisee);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_binarisee;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Inversion_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Inverser les niveaux de gris de l'image
    image_inversee = 255 - image;

    % Convertir l'image inversée en uint8 pour l'affichage
    image_inversee = uint8(image_inversee);

    % Afficher l'image inversée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_inversee);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_inversee;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Luminosite_Callback(hObject, eventdata, handles)
% hObject    handle to Luminosite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Contraste_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Déterminer les valeurs minimales et maximales de l'image
    min_val = min(image(:));
    max_val = max(image(:));

    % Appliquer une transformation pour augmenter le contraste
    image_contraste_plus = 255 * (image - min_val) / (max_val - min_val);

    % Convertir l'image en uint8 pour l'affichage
    image_contraste_plus = uint8(image_contraste_plus);

    % Afficher l'image avec un contraste accru dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_contraste_plus);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_contraste_plus;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function EgalisationDeLhistogramme_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Égaliser l'histogramme de l'image
    image_egalisee = histeq(uint8(image));

    % Afficher l'image égalisée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_egalisee);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_egalisee;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Vérifier si l'image est en couleur
    if size(image, 3) == 3
        % L'image est en couleur, convertir en niveaux de gris
        image_gray = rgb2gray(image);

        % Calcul de l'histogramme pour chaque canal de couleur
        [counts_red, bins_red] = imhist(image(:,:,1));
        [counts_green, bins_green] = imhist(image(:,:,2));
        [counts_blue, bins_blue] = imhist(image(:,:,3));

        % Tracer l'histogramme
        axes(handles.imgS);
        hold off; % Pour effacer les anciens graphiques s'il y en a
        bar(bins_red, counts_red, 'BarWidth', 0.8, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        hold on;
        bar(bins_green, counts_green, 'BarWidth', 0.8, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        bar(bins_blue, counts_blue, 'BarWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        % Mettre en place les étiquettes
        xlabel('Niveau de gris');
        ylabel('Nombre de pixels');
        title('Histogramme de l''image (Couleurs)');

        % Réglage de l'axe x pour couvrir la plage complète des niveaux de gris
        xlim([0, 255]);

        % Affichage de la grille
        grid on;
        legend('Rouge', 'Vert', 'Bleu');
        hold off;
    else
        % L'image est en niveaux de gris, calculer l'histogramme directement
        [counts, bins] = imhist(image);

        % Tracer l'histogramme
        axes(handles.imgS);
        bar(bins, counts, 'BarWidth', 1);

        % Mettre en place les étiquettes
        xlabel('Niveau de gris');
        ylabel('Nombre de pixels');
        title('Histogramme de l''image (Niveaux de gris)');

        % Réglage de l'axe x pour couvrir la plage complète des niveaux de gris
        xlim([0, 255]);

        % Affichage de la grille
        grid on;
    end


% --------------------------------------------------------------------
function SUSAN_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % =======================conversion de l'image=============================
    d = length(size(image));
    if d == 3
        image_gray = double(rgb2gray(image));
    elseif d == 2
        image_gray = double(image);
    end

    [n, m] = size(image_gray);
    % =============================données=====================================
    rayon = 1;
    alpha = 80;
    r = 5;
    alpha = alpha / 100;

    % ========================génerateur de mask=============================
    mask = zeros(2 * rayon + 1);
    b = ones(rayon + 1);
    for i = 1:rayon + 1
        for j = 1:rayon + 1
            if (rayon == 1)
                if(j > i)
                    b(i, j) = 0;
                end
            else
                if(j > i + 1)
                    b(i, j) = 0;
                end
            end
        end
    end
    mask(1:rayon + 1, rayon + 1:2 * rayon + 1) = b;
    mask(1:rayon + 1, 1:rayon + 1) = rot90(b);
    mask0 = mask;
    mask0 = flipdim(mask0, 1);
    mask = mask0 + mask;
    mask(rayon + 1, :) = mask(rayon + 1, :) - 1;

    % ==========================réponse maximale============================
    max_reponse = sum(sum(mask));
    % =====================balayage de toute l'image===========================
    f = zeros(n, m);
    for i = (rayon + 1):n - rayon
        for j = (rayon + 1):m - rayon
            image_courant = image_gray(i - rayon:i + rayon, j - rayon:j + rayon);
            image_courant_mask = image_courant .* mask;
            inteniste_cental = image_courant_mask(rayon + 1, rayon + 1);
            s = exp(-1 * (((image_courant_mask - inteniste_cental) / max_reponse) .^ 6));
            somme = sum(sum(s));
            % si le centre du mask est un 0 il faut soustraire les zeros des filtres
            if (inteniste_cental == 0)
                somme = somme - length((find(mask == 0)));
            end
            f(i, j) = somme;
        end
    end

    % =============selection et seuillage des points d'interét=================
    ff = f(rayon + 1:n - (rayon + 1), rayon + 1:m - (rayon + 1));
    minf = min(min(ff));
    maxf = max(max(f));
    fff = f;
    d = 2 * r + 1;
    temp1 = round(n / d);
    if (temp1 - n / d) < 0.5 && (temp1 - n / d) > 0
        temp1 = temp1 - 1;
    end
    temp2 = round(m / d);
    if (temp2 - m / d) < 0.5 && (temp2 - m / d) > 0
        temp2 = temp2 - 1;
    end
    fff(n:temp1 * d + d, m:temp2 * d + d) = 0;

    for i = (r + 1):d:temp1 * d + d
        for j = (r + 1):d:temp2 * d + d
            window = fff(i - r:i + r, j - r:j + r);
            window0 = window;
            [xx, yy] = find(window0 == 0);
            for k = 1:length(xx)
                window0(xx(k), yy(k)) = max(max(window0));
            end
            minwindow = min(min(window0));
            [y, x] = find(minwindow ~= window & window <= minf + alpha * (maxf - minf) & window > 0);
            [u, v] = find(minwindow == window);
            if length(u) > 1
                for l = 2:length(u)
                    fff(i - r - 1 + u(l), j - r - 1 + v(l)) = 0 ;
                end
            end
            if length(x) ~= 0
                for l = 1:length(y)
                    fff(i - r - 1 + y(l), j - r - 1 + x(l)) = 0 ;
                end
            end
        end
    end
    seuil = minf + alpha * (maxf - minf);
    [u, v] = find(minf <= fff & fff <= seuil );

    % ==============affichage des resultats====================================
    axes(handles.imgS);
    imshow(image);
    hold on;
    plot(v, u, '.r', 'MarkerSize', 10);
    hold off;
    message = sprintf(' le nombre des points d''intérêts: %d      ', length(v));
    msgbox(message);



% --------------------------------------------------------------------
function HARRIS_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Conversion de l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image_gray = rgb2gray(image);
    else
        image_gray = image;
    end

    % Paramètres pour la détection de coins de Harris
    lambda = 0.04;
    sigma = 1;
    seuil = 200;
    r = 6;
    w = 5 * sigma;

    % Calcul des dérivées horizontales et verticales
    dx = [-1 0 1;
          -2 0 2;
          -1 0 1]; % Filtre de Sobel pour la dérivée horizontale
    dy = dx'; % Transposée pour la dérivée verticale
    g = fspecial('gaussian', max(1, fix(w)), sigma); % Filtre Gaussien

    % Calcul des dérivées selon x et y
    Ix = conv2(double(image_gray), dx, 'same');
    Iy = conv2(double(image_gray), dy, 'same');

    % Calcul des composantes de la matrice de Harris
    Ix2 = conv2(Ix .^ 2, g, 'same');
    Iy2 = conv2(Iy .^ 2, g, 'same');
    Ixy = conv2(Ix .* Iy, g, 'same');
    detM = Ix2 .* Iy2 - Ixy .^ 2;
    trM = Ix2 + Iy2;
    R = detM - lambda * trM .^ 2;

    % Seuillage et suppression des non-maxima locaux
    R1 = (1000 / (1 + max(max(R)))) * R;
    [u, v] = find(R1 <= seuil);
    nb = length(u);
    for k = 1:nb
        R1(u(k), v(k)) = 0;
    end

    % Suppression des non-maxima locaux
    [m, n] = size(R1);
    R11 = zeros(m + 2 * r, n + 2 * r);
    R11(r + 1:m + r, r + 1:n + r) = R1;
    [m1, n1] = size(R11);
    for i = r + 1:m1 - r
        for j = r + 1:n1 - r
            fenetre = R11(i - r:i + r, j - r:j + r);
            ma = max(max(fenetre));
            if fenetre(r + 1, r + 1) < ma
                R11(i, j) = 0;
            end
        end
    end
    R11 = R11(r + 1:m + r, r + 1:n + r);
    [x, y] = find(R11);

    % Affichage des résultats
    axes(handles.imgS);
    imshow(image);
    hold on;
    plot(y, x, '.r', 'MarkerSize', 10);
    hold off;
    message = sprintf('Le nombre des points d''intérêts détectés : %d', length(x));
    msgbox(message);


% --------------------------------------------------------------------
function HARRIS_ModeleElectrostatique_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Conversion de l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image_gray = rgb2gray(image);
    else
        image_gray = image;
    end

    % Paramètres pour la détection de points d'intérêt
    k = 0.04;
    sigma = 1;
    seuil = 100;
    r = 6;

    % Conversion en double pour les calculs
    imd = double(image_gray);

    % Calcul des dérivées horizontales et verticales
    dxa = [-sqrt(2)/4, 0, sqrt(2)/4; -1, 0, 1; -sqrt(2)/4, 0, sqrt(2)/4];
    dya = dxa'; % Transposée pour la dérivée verticale

    % Filtre gaussien
    g = fspecial('gaussian', max(1, fix(5 * sigma)), sigma);

    % Calcul des dérivées selon x et y
    Ixa = conv2(imd, dxa, 'same');
    Iya = conv2(imd, dya, 'same');

    % Calcul des produits de dérivées
    Ixa2 = conv2(Ixa .^ 2, g, 'same');
    Iya2 = conv2(Iya .^ 2, g, 'same');
    Ixya = conv2(Ixa .* Iya, g, 'same');

    % Calcul de la réponse de Harris
    R = Ixa2 .* Iya2 - Ixya .^ 2 - k * (Ixa2 + Iya2) .^ 2;

    % Normalisation de la réponse
    R1 = (1000 / (max(max(R)))) * R;

    % Seuillage et suppression des non-maxima locaux
    [u, v] = find(R1 <= seuil);
    nb = length(u);
    for k = 1:nb
        R1(u(k), v(k)) = 0;
    end

    % Suppression des non-maxima locaux
    [m, n] = size(R1);
    R11 = zeros(m + 2 * r, n + 2 * r);
    R11(r + 1:m + r, r + 1:n + r) = R1;
    [m1, n1] = size(R11);
    for i = r + 1:m1 - r
        for j = r + 1:n1 - r
            fenetre = R11(i - r:i + r, j - r:j + r);
            ma = max(max(fenetre));
            if fenetre(r + 1, r + 1) < ma
                R11(i, j) = 0;
            end
        end
    end
    R11 = R11(r + 1:m + r, r + 1:n + r);
    [x, y] = find(R11);

    % Affichage des résultats
    axes(handles.imgS);
    imshow(image);
    hold on;
    plot(y, x, '.r', 'MarkerSize', 10);
    hold off;
    message = sprintf('Le nombre des points d''intérêts détectés : %d', length(x));
    msgbox(message);


% --------------------------------------------------------------------
function Erosion_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour l'érosion
    se = ones(3); % Élément structurant 3x3

    % Taille de l'image
    [m, n] = size(image);

    % Initialiser une nouvelle image pour stocker le résultat de l'érosion
    image_eroded = zeros(size(image));

    % Appliquer l'érosion manuellement
    for i = 2:m-1
        for j = 2:n-1
            % Extraire la région de l'image sous l'élément structurant
            region = image(i-1:i+1, j-1:j+1);
            
            % Calculer le minimum de la convolution
            min_value = min(region(se > 0));
            
            % Assigner la valeur minimale à l'image érodée
            image_eroded(i, j) = min_value;
        end
    end

    % Afficher l'image érodée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(image_eroded));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_eroded;
    guidata(hObject, handles);

% --------------------------------------------------------------------
function Dilatation_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour la dilatation
    se = ones(3); % Élément structurant 3x3

    % Taille de l'image
    [m, n] = size(image);

    % Initialiser une nouvelle image pour stocker le résultat de la dilatation
    image_dilated = zeros(size(image));

    % Appliquer la dilatation manuellement
    for i = 2:m-1
        for j = 2:n-1
            % Extraire la région de l'image sous l'élément structurant
            region = image(i-1:i+1, j-1:j+1);
            
            % Calculer le maximum de la convolution
            max_value = max(region(se > 0));
            
            % Assigner la valeur maximale à l'image dilatée
            image_dilated(i, j) = max_value;
        end
    end

    % Afficher l'image dilatée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(image_dilated));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_dilated;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Ouverture_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour l'ouverture
    se = strel('disk', 3); % Élément structurant de forme circulaire avec un rayon de 3 pixels

    % Appliquer l'érosion à l'image
    image_erodee = imerode(image, se);

    % Appliquer la dilatation à l'image érodée
    image_ouverte = imdilate(image_erodee, se);

    % Afficher l'image ouverte dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_ouverte);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_ouverte;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Fermeture_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour la fermeture
    se = strel('disk', 3); % Élément structurant de forme circulaire avec un rayon de 3 pixels

    % Appliquer la dilatation à l'image
    image_dilatee = imdilate(image, se);

    % Appliquer l'érosion à l'image dilatée
    image_fermee = imerode(image_dilatee, se);

    % Afficher l'image fermée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_fermee);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_fermee;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Contour_Callback(hObject, eventdata, handles)
% hObject    handle to Contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Spectre_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs
    
    % Calculer la transformée de Fourier 2D (FFT) de l'image
    spectrum = fft2(image);

    % Décaler le spectre pour que les basses fréquences soient au centre
    spectrum_shifted = fftshift(spectrum);

    % Calculer le module du spectre (amplitude)
    spectrum_magnitude = abs(spectrum_shifted);

    % Échelle logarithmique pour une meilleure visualisation
    spectrum_log = log(1 + spectrum_magnitude);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    spectrum_normalized = uint8(255 * mat2gray(spectrum_log));

    % Afficher l'image du spectre dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(spectrum_normalized, []);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = spectrum_normalized;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function FiltragePasseBasIdeal_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Séparation des canaux de couleur si l'image est en couleur
    if size(image, 3) == 3
        % Canal rouge
        image_r = image(:,:,1);
        % Canal vert
        image_g = image(:,:,2);
        % Canal bleu
        image_b = image(:,:,3);

        % Filtrage pour chaque canal de couleur
        filtered_image_r = apply_filterBasIdeal(image_r);
        filtered_image_g = apply_filterBasIdeal(image_g);
        filtered_image_b = apply_filterBasIdeal(image_b);

        % Reconstruction de l'image filtrée
        filtered_image = cat(3, filtered_image_r, filtered_image_g, filtered_image_b);
    else
        % Filtrage direct si l'image est en niveaux de gris
        filtered_image = apply_filterBasIdeal(image);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(filtered_image));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);

function filtered_image = apply_filterBasIdeal(image)
    % Calculer la transformée de Fourier 2D de l'image
    F = fft2(image);

    % Décaler le spectre pour que les basses fréquences soient au centre
    F_shifted = fftshift(F);

    % Taille de l'image
    [M, N] = size(image);
    % Coordonnées du centre de l'image
    u0 = ceil(M/2);
    v0 = ceil(N/2);

    % Définir le rayon du filtre passe-bas idéal
    D0 = 30; % Vous pouvez ajuster cette valeur selon vos besoins

    % Créer le masque du filtre passe-bas idéal
    H = zeros(M, N);
    for u = 1:M
        for v = 1:N
            D_uv = sqrt((u - u0)^2 + (v - v0)^2);
            if D_uv <= D0
                H(u, v) = 1;
            end
        end
    end

    % Appliquer le filtre en multipliant le spectre par le masque
    G_shifted = F_shifted .* H;

    % Décaler inversement le spectre
    G = ifftshift(G_shifted);

    % Calculer l'image filtrée en prenant la transformée de Fourier inverse
    filtered_image = abs(ifft2(G));

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    filtered_image = uint8(255 * mat2gray(filtered_image));

% --------------------------------------------------------------------
function FiltragePasseHautIdeal_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Séparation des canaux de couleur si l'image est en couleur
    if size(image, 3) == 3
        % Canal rouge
        image_r = image(:,:,1);
        % Canal vert
        image_g = image(:,:,2);
        % Canal bleu
        image_b = image(:,:,3);

        % Filtrage pour chaque canal de couleur
        filtered_image_r = apply_filterHautIdeal(image_r);
        filtered_image_g = apply_filterHautIdeal(image_g);
        filtered_image_b = apply_filterHautIdeal(image_b);

        % Reconstruction de l'image filtrée
        filtered_image = cat(3, filtered_image_r, filtered_image_g, filtered_image_b);
    else
        % Filtrage direct si l'image est en niveaux de gris
        filtered_image = apply_filterHautIdeal(image);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(filtered_image));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


function filtered_image = apply_filterHautIdeal(image)
    % Calculer la transformée de Fourier 2D de l'image
    F = fft2(image);

    % Décaler le spectre pour que les basses fréquences soient au centre
    F_shifted = fftshift(F);

    % Taille de l'image
    [M, N] = size(image);
    % Coordonnées du centre de l'image
    u0 = ceil(M/2);
    v0 = ceil(N/2);

    % Définir le rayon du filtre passe-haut idéal
    D0 = 30; % Vous pouvez ajuster cette valeur selon vos besoins

    % Créer le masque du filtre passe-haut idéal
    H = ones(M, N);
    for u = 1:M
        for v = 1:N
            D_uv = sqrt((u - u0)^2 + (v - v0)^2);
            if D_uv <= D0
                H(u, v) = 0;
            end
        end
    end

    % Appliquer le filtre en multipliant le spectre par le masque
    G_shifted = F_shifted .* H;

    % Décaler inversement le spectre
    G = ifftshift(G_shifted);

    % Calculer l'image filtrée en prenant la transformée de Fourier inverse
    filtered_image = abs(ifft2(G));

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    filtered_image = uint8(255 * mat2gray(filtered_image));


% --------------------------------------------------------------------
function FiltragePasseBasBUTTERWORTH_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Séparation des canaux de couleur si l'image est en couleur
    if size(image, 3) == 3
        % Canal rouge
        image_r = image(:,:,1);
        % Canal vert
        image_g = image(:,:,2);
        % Canal bleu
        image_b = image(:,:,3);

        % Filtrage pour chaque canal de couleur
        filtered_image_r = apply_filterBasBUTTERWORTH(image_r);
        filtered_image_g = apply_filterBasBUTTERWORTH(image_g);
        filtered_image_b = apply_filterBasBUTTERWORTH(image_b);

        % Reconstruction de l'image filtrée
        filtered_image = cat(3, filtered_image_r, filtered_image_g, filtered_image_b);
    else
        % Filtrage direct si l'image est en niveaux de gris
        filtered_image = apply_filterBasBUTTERWORTH(image);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(filtered_image));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


function filtered_image = apply_filterBasBUTTERWORTH(image)
    % Calculer la transformée de Fourier 2D de l'image
    F = fft2(image);

    % Décaler le spectre pour que les basses fréquences soient au centre
    F_shifted = fftshift(F);

    % Taille de l'image
    [M, N] = size(image);
    % Coordonnées du centre de l'image
    u0 = ceil(M/2);
    v0 = ceil(N/2);

    % Définir la fréquence de coupure et l'ordre du filtre Butterworth
    D0 = 30; % Fréquence de coupure (ajustez selon vos besoins)
    n = 2;   % Ordre du filtre Butterworth (ajustez selon vos besoins)

    % Créer le masque du filtre passe-bas Butterworth
    H = zeros(M, N);
    for u = 1:M
        for v = 1:N
            D_uv = sqrt((u - u0)^2 + (v - v0)^2);
            H(u, v) = 1 / (1 + (D_uv / D0)^(2*n));
        end
    end

    % Appliquer le filtre en multipliant le spectre par le masque
    G_shifted = F_shifted .* H;

    % Décaler inversement le spectre
    G = ifftshift(G_shifted);

    % Calculer l'image filtrée en prenant la transformée de Fourier inverse
    filtered_image = abs(ifft2(G));

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    filtered_image = uint8(255 * mat2gray(filtered_image));


% --------------------------------------------------------------------
function FiltragePasseHautBUTTERWORTH_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Séparation des canaux de couleur si l'image est en couleur
    if size(image, 3) == 3
        % Canal rouge
        image_r = image(:,:,1);
        % Canal vert
        image_g = image(:,:,2);
        % Canal bleu
        image_b = image(:,:,3);

        % Filtrage pour chaque canal de couleur
        filtered_image_r = apply_filterHautBUTTERWORTH(image_r);
        filtered_image_g = apply_filterHautBUTTERWORTH(image_g);
        filtered_image_b = apply_filterHautBUTTERWORTH(image_b);

        % Reconstruction de l'image filtrée
        filtered_image = cat(3, filtered_image_r, filtered_image_g, filtered_image_b);
    else
        % Filtrage direct si l'image est en niveaux de gris
        filtered_image = apply_filterHautBUTTERWORTH(image);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(uint8(filtered_image));

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


function filtered_image = apply_filterHautBUTTERWORTH(image)
    % Calculer la transformée de Fourier 2D de l'image
    F = fft2(image);

    % Décaler le spectre pour que les basses fréquences soient au centre
    F_shifted = fftshift(F);

    % Taille de l'image
    [M, N] = size(image);
    % Coordonnées du centre de l'image
    u0 = ceil(M/2);
    v0 = ceil(N/2);

    % Définir la fréquence de coupure et l'ordre du filtre Butterworth
    D0 = 30; % Fréquence de coupure (ajustez selon vos besoins)
    n = 2;   % Ordre du filtre Butterworth (ajustez selon vos besoins)

    % Créer le masque du filtre passe-haut Butterworth
    H = ones(M, N);
    for u = 1:M
        for v = 1:N
            D_uv = sqrt((u - u0)^2 + (v - v0)^2);
            H(u, v) = 1 / (1 + (D0 / D_uv)^(2*n));
        end
    end

    % Appliquer le filtre en multipliant le spectre par le masque
    G_shifted = F_shifted .* H;

    % Décaler inversement le spectre
    G = ifftshift(G_shifted);

    % Calculer l'image filtrée en prenant la transformée de Fourier inverse
    filtered_image = abs(ifft2(G));

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    filtered_image = uint8(255 * mat2gray(filtered_image));



% --------------------------------------------------------------------
function FiltrePasseBas_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseBas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltrePasseHaut_Callback(hObject, eventdata, handles)
% hObject    handle to FiltrePasseHaut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltreGradient_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir les opérateurs de Sobel pour le gradient
    sobel_x = [-1 0 1; -2 0 2; -1 0 1];
    sobel_y = [-1 -2 -1; 0 0 0; 1 2 1];

    % Appliquer les filtres de Sobel pour obtenir les gradients en x et y pour chaque canal
    [n, m, d] = size(image);
    gradient_x = zeros(n, m, d);
    gradient_y = zeros(n, m, d);

    for channel = 1:d
        gradient_x(:, :, channel) = conv2(image(:, :, channel), sobel_x, 'same');
        gradient_y(:, :, channel) = conv2(image(:, :, channel), sobel_y, 'same');
    end

    % Calculer le gradient total (magnitude) pour chaque canal
    gradient_magnitude = sqrt(gradient_x.^2 + gradient_y.^2);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255 pour chaque canal
    gradient_magnitude = uint8(255 * mat2gray(gradient_magnitude));

    % Afficher l'image du gradient dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(gradient_magnitude);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = gradient_magnitude;
    guidata(hObject, handles);

% --------------------------------------------------------------------
function FiltreLaplacien_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir le noyau du filtre Laplacien
    laplacian_kernel = [0, -1, 0; -1, 4, -1; 0, -1, 0];

    % Appliquer le filtre laplacien à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(laplacian_kernel, image(:, :, channel), 'same');
    end

    % Normaliser les valeurs pour obtenir une image entre 0 et 255 pour chaque canal
    filtered_image = uint8(255 * mat2gray(filtered_image));

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreSobel_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir les opérateurs de Sobel pour le gradient en x et y
    sobel_x = [-1 0 1; -2 0 2; -1 0 1];
    sobel_y = [-1 -2 -1; 0 0 0; 1 2 1];

    % Appliquer les filtres de Sobel pour obtenir les gradients en x et y pour chaque canal
    [n, m, d] = size(image);
    gradient_x = zeros(n, m, d);
    gradient_y = zeros(n, m, d);

    for channel = 1:d
        gradient_x(:, :, channel) = conv2(image(:, :, channel), sobel_x, 'same');
        gradient_y(:, :, channel) = conv2(image(:, :, channel), sobel_y, 'same');
    end

    % Calculer la magnitude du gradient pour chaque canal
    gradient_magnitude = abs(gradient_x) + abs(gradient_y);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255 pour chaque canal
    gradient_magnitude = uint8(255 * mat2gray(gradient_magnitude));

    % Afficher l'image du gradient dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(gradient_magnitude);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = gradient_magnitude;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltrePrewitt_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir les opérateurs de Prewitt pour le gradient en x et y
    prewitt_x = [-1 0 1; -1 0 1; -1 0 1];
    prewitt_y = [-1 -1 -1; 0 0 0; 1 1 1];

    % Appliquer les filtres de Prewitt pour obtenir les gradients en x et y pour chaque canal
    [n, m, d] = size(image);
    gradient_x = zeros(n, m, d);
    gradient_y = zeros(n, m, d);

    for channel = 1:d
        gradient_x(:, :, channel) = conv2(image(:, :, channel), prewitt_x, 'same');
        gradient_y(:, :, channel) = conv2(image(:, :, channel), prewitt_y, 'same');
    end

    % Calculer la magnitude du gradient pour chaque canal
    gradient_magnitude = abs(gradient_x) + abs(gradient_y);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255 pour chaque canal
    gradient_magnitude = uint8(255 * mat2gray(gradient_magnitude));

    % Afficher l'image du gradient dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(gradient_magnitude);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = gradient_magnitude;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function FiltreRobet_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir les opérateurs de Roberts pour le gradient en x et y
    roberts_x = [1 0; 0 -1];
    roberts_y = [0 1; -1 0];

    % Appliquer les filtres de Roberts pour obtenir les gradients en x et y pour chaque canal
    [n, m, d] = size(image);
    gradient_x = zeros(n, m, d);
    gradient_y = zeros(n, m, d);

    for channel = 1:d
        gradient_x(:, :, channel) = conv2(image(:, :, channel), roberts_x, 'same');
        gradient_y(:, :, channel) = conv2(image(:, :, channel), roberts_y, 'same');
    end

    % Calculer la magnitude du gradient pour chaque canal
    gradient_magnitude = abs(gradient_x) + abs(gradient_y);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255 pour chaque canal
    gradient_magnitude = uint8(255 * mat2gray(gradient_magnitude));

    % Afficher l'image du gradient dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(gradient_magnitude);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = gradient_magnitude;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function FiltreKirsch_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir les 8 masques de Kirsch
    kirsch_masks = {
        [-3 -3 5; -3 0 5; -3 -3 5],
        [-3 5 5; -3 0 5; -3 -3 -3],
        [5 5 5; -3 0 -3; -3 -3 -3],
        [5 5 -3; 5 0 -3; -3 -3 -3],
        [5 -3 -3; 5 0 -3; 5 -3 -3],
        [-3 -3 -3; 5 0 -3; 5 5 -3],
        [-3 -3 -3; -3 0 -3; 5 5 5],
        [-3 -3 -3; -3 0 5; -3 5 5]
    };

    % Appliquer les masques de Kirsch et trouver le maximum des magnitudes de gradient
    gradient_magnitude = zeros(size(image));
    for i = 1:length(kirsch_masks)
        mask = kirsch_masks{i};
        gradient = zeros(size(image, 1), size(image, 2), size(image, 3));
        for channel = 1:size(image, 3)
            gradient(:, :, channel) = conv2(image(:, :, channel), mask, 'same');
        end
        gradient_magnitude = max(gradient_magnitude, abs(sum(gradient, 3)));
    end

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    gradient_magnitude = uint8(255 * mat2gray(gradient_magnitude));

    % Afficher l'image du gradient dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(gradient_magnitude);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = gradient_magnitude;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreMarrHildreth_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs
    
    % Paramètres du filtre Marr-Hildreth
    sigma = 1.5; % Écart type du filtre gaussien

    % Appliquer le lissage gaussien
    smoothed_image = imgaussfilt(image, sigma);

    % Appliquer l'opérateur laplacien pour obtenir le LoG
    laplacian_of_gaussian = del2(smoothed_image);

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    LoG_normalized = uint8(255 * mat2gray(laplacian_of_gaussian));

    % Afficher l'image du LoG dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(LoG_normalized);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = LoG_normalized;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Lineaires_Callback(hObject, eventdata, handles)
% hObject    handle to Lineaires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NonLineaires_Callback(hObject, eventdata, handles)
% hObject    handle to NonLineaires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ContourInterieur_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour l'érosion
    se = strel('disk', 3); % Élément structurant de forme circulaire avec un rayon de 3 pixels

    % Appliquer une érosion à l'image
    image_erodee = imerode(image, se);

    % Soustraire l'image érodée de l'image originale pour obtenir le contour intérieur
    contour_interieur = image - image_erodee;

    % Afficher le contour intérieur dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(contour_interieur);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = contour_interieur;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function ContourExterieur_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour la dilatation
    se = strel('disk', 3); % Élément structurant de forme circulaire avec un rayon de 3 pixels

    % Appliquer une dilatation à l'image
    image_dilatee = imdilate(image, se);

    % Soustraire l'image originale de l'image dilatée pour obtenir le contour extérieur
    contour_exterieur = image_dilatee - image;

    % Afficher le contour extérieur dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(contour_exterieur);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = contour_exterieur;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function ContourMorphologique_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Définir l'élément structurant pour l'opération morphologique
    se = strel('disk', 3); % Élément structurant de forme circulaire avec un rayon de 3 pixels

    % Appliquer une dilatation suivie d'une érosion à l'image
    image_dilatee = imdilate(image, se);
    image_erodee = imerode(image, se);

    % Soustraire l'image érodée de l'image dilatée pour obtenir le contour morphologique
    contour_morphologique = image_dilatee - image_erodee;

    % Afficher le contour morphologique dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(contour_morphologique);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = contour_morphologique;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function LPlus_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Augmenter la luminosité en ajoutant une valeur constante (ajustez selon vos besoins)
    valeur_constante = 50;
    image_luminosite_plus = image + valeur_constante;

    % Assurer que les valeurs de l'image restent dans la plage [0, 255]
    image_luminosite_plus(image_luminosite_plus > 255) = 255;

    % Convertir l'image en uint8 pour l'affichage
    image_luminosite_plus = uint8(image_luminosite_plus);

    % Afficher l'image avec une luminosité augmentée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_luminosite_plus);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_luminosite_plus;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function LMoins_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Diminuer la luminosité en soustrayant une valeur constante (ajustez selon vos besoins)
    valeur_constante = 50;
    image_luminosite_moins = image - valeur_constante;

    % Assurer que les valeurs de l'image restent dans la plage [0, 255]
    image_luminosite_moins(image_luminosite_moins < 0) = 0;

    % Convertir l'image en uint8 pour l'affichage
    image_luminosite_moins = uint8(image_luminosite_moins);

    % Afficher l'image avec une luminosité diminuée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_luminosite_moins);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_luminosite_moins;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreMedian_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function FiltreMoyenneur_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreMoyenneur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltreGaussien_Callback(hObject, eventdata, handles)
% hObject    handle to FiltreGaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FiltrePyramidal5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir le noyau du filtre pyramidal 5x5
    kernel = (1/81) * [1, 2, 3, 2, 1; 2, 4, 6, 4, 2; 3, 6, 9, 6, 3; 2, 4, 6, 4, 2; 1, 2, 3, 2, 1];

    % Appliquer le filtre pyramidal à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(kernel, image(:, :, channel), 'same');
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);




% --------------------------------------------------------------------
function FiltreConique5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir le noyau du filtre conique 5x5
    kernel = (1/25) * [0, 0, 1, 0, 0; 0, 2, 2, 2, 0; 1, 2, 5, 2, 1; 0, 2, 2, 2, 0; 0, 0, 1, 0, 0];

    % Appliquer le filtre conique à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(kernel, image(:, :, channel), 'same');
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreBinomial5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Définir le noyau du filtre binomial 5x5
    kernel = (1/256) * [1, 4, 6, 4, 1; 4, 16, 24, 16, 4; 6, 24, 36, 24, 6; 4, 16, 24, 16, 4; 1, 4, 6, 4, 1];

    % Appliquer le filtre binomial à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(kernel, image(:, :, channel), 'same');
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreGaussien3x3_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Créer le noyau de filtre gaussien 3x3
    kernel = 1/16 * [1, 2, 1; 2, 4, 2; 1, 2, 1];

    % Appliquer le filtre gaussien à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(kernel, image(:, :, channel), 'same');
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreGaussien5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Créer le noyau de filtre gaussien 5x5
    kernel = 1/273 * [1, 4, 7, 4, 1; 4, 16, 26, 16, 4; 7, 26, 41, 26, 7; 4, 16, 26, 16, 4; 1, 4, 7, 4, 1];

    % Appliquer le filtre gaussien à chaque canal de l'image
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        filtered_image(:, :, channel) = filter2(kernel, image(:, :, channel), 'same');
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreMoyenneur3x3_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Appliquer le filtre moyenneur 3x3
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        for x = 2 : n-1
            for y = 2 : m-1
                subimage = image(x-1:x+1, y-1:y+1, channel);
                filtered_image(x, y, channel) = mean(subimage(:));
            end
        end
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreMoyenneur5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Appliquer le filtre moyenneur 5x5
    [n, m, d] = size(image);
    filtered_image = zeros(n, m, d);

    for channel = 1:d
        for x = 3 : n-2
            for y = 3 : m-2
                subimage = image(x-2:x+2, y-2:y+2, channel);
                filtered_image(x, y, channel) = mean(subimage(:));
            end
        end
    end

    % Convertir l'image filtrée en type uint8 pour l'affichage
    filtered_image = uint8(filtered_image);

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function FiltreMedian3x3_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data; % Supposons que votre image soit stockée dans handles.courant_data

    % Vérifier si l'image est en niveaux de gris ou en couleur
    if size(image, 3) == 3
        % L'image est en couleur
        % Appliquer le filtre médian 3x3 sur chaque canal de l'image
        filtered_image(:,:,1) = medfilt2(image(:,:,1), [3, 3]);
        filtered_image(:,:,2) = medfilt2(image(:,:,2), [3, 3]);
        filtered_image(:,:,3) = medfilt2(image(:,:,3), [3, 3]);
    else
        % L'image est en niveaux de gris
        % Appliquer le filtre médian 3x3
        filtered_image = medfilt2(image, [3, 3]);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function FiltreMedian5x5_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data; % Supposons que votre image soit stockée dans handles.courant_data

    % Vérifier si l'image est en niveaux de gris ou en couleur
    if size(image, 3) == 3
        % L'image est en couleur
        % Appliquer le filtre médian 5x5 sur chaque canal de l'image
        filtered_image(:,:,1) = medfilt2(image(:,:,1), [5, 5]);
        filtered_image(:,:,2) = medfilt2(image(:,:,2), [5, 5]);
        filtered_image(:,:,3) = medfilt2(image(:,:,3), [5, 5]);
    else
        % L'image est en niveaux de gris
        % Appliquer le filtre médian 5x5
        filtered_image = medfilt2(image, [5, 5]);
    end

    % Afficher l'image filtrée dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(filtered_image);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = filtered_image;
    guidata(hObject, handles);



% --------------------------------------------------------------------
function Bruit_Callback(hObject, eventdata, handles)
% hObject    handle to Bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gaussien_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Taille de l'image
    [M, N, ~] = size(image);

    % Paramètres du bruit gaussien
    mu = 0;    % Moyenne
    sigma = 20; % Écart type (ajustez selon vos besoins)

    % Générer un bruit gaussien pour chaque canal
    noise = sigma * randn(M, N, 3) + mu;

    % Ajouter le bruit gaussien à l'image
    image_noise = image + noise;

    % Normaliser les valeurs pour obtenir une image entre 0 et 255
    image_noise = uint8(255 * mat2gray(image_noise));

    % Afficher l'image avec bruit gaussien dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_noise);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_noise;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Proivre_et_sel_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = double(handles.courant_data); % Convertir en double pour les calculs

    % Taille de l'image
    [M, N] = size(image);

    % Paramètres du bruit de poivre et sel
    density = 0.05; % Densité du bruit (ajustez selon vos besoins)

    % Générer un masque aléatoire pour le bruit de poivre et sel
    mask = rand(M, N);

    % Appliquer le bruit de poivre et sel à l'image
    image_salt_pepper = image;
    image_salt_pepper(mask <= density/2) = 0;     % Poivre
    image_salt_pepper(mask > 1 - density/2) = 255; % Sel

    % Convertir l'image en uint8 pour l'affichage
    image_salt_pepper = uint8(image_salt_pepper);

    % Afficher l'image avec bruit de poivre et sel dans l'axe spécifié dans handles (handles.imgS)
    axes(handles.imgS);
    imshow(image_salt_pepper);

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image_salt_pepper;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function Hough_Callback(hObject, eventdata, handles)
% hObject    handle to Hough (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HoughDroites_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Conversion de l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image_gray = rgb2gray(image);
    else
        image_gray = image;
    end

    % Détection de contours avec l'opérateur de Canny
    edges = edge(image_gray, 'Canny');

    % Calcul de la transformée de Hough
    [H, theta, rho] = hough(edges);

    % Recherche des pics dans la transformée de Hough
    P = houghpeaks(H, 5, 'threshold', ceil(0.3 * max(H(:))));

    % Extraction des coordonnées des droites
    lines = houghlines(edges, theta, rho, P);

    % Affichage des résultats
    axes(handles.imgS);
    imshow(image);
    hold on;

    % Affichage des droites détectées
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'r');
        plot(xy(1, 1), xy(1, 2), 'x', 'LineWidth', 2, 'Color', 'g');
        plot(xy(2, 1), xy(2, 2), 'x', 'LineWidth', 2, 'Color', 'g');
    end

    hold off;

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image;
    guidata(hObject, handles);


% --------------------------------------------------------------------
function HoughCercles_Callback(hObject, eventdata, handles)
    % Obtenir l'image à partir des données stockées dans handles
    image = handles.courant_data;

    % Conversion de l'image en niveaux de gris si elle est en couleur
    if size(image, 3) == 3
        image_gray = rgb2gray(image);
    else
        image_gray = image;
    end

    % Réduction du bruit avec un filtre gaussien
    image_blurred = imgaussfilt(image_gray, 2);

    % Détection de bords avec l'opérateur de Canny
    edges = edge(image_blurred, 'Canny');

    % Paramètres de la transformée de Hough pour les cercles
    rayonMin = 10;
    rayonMax = 80;

    % Détection des cercles avec la transformée de Hough
    [centers, radii, ~] = imfindcircles(edges, [rayonMin, rayonMax]);

    % Affichage des résultats
    axes(handles.imgS);
    imshow(image);
    hold on;

    % Affichage des cercles détectés
    viscircles(centers, radii, 'EdgeColor', 'b');

    hold off;

    % Mettre à jour les données de l'interface utilisateur si nécessaire
    handles.ima_traite = image;
    guidata(hObject, handles);
