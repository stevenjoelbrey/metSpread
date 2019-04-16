class relaxo() :
    """
    Custom implementation of the double lasso. The idea is to use lasso twice. 
    First time through use it explicitly for variable selection (what features
    have coef values not equal to zero) and then again for linear regression
    fitting regularization. The second use only uses features that had non-zero
    coefs on X% of external loop leave one out iterations. 

    """
    
    def __init__(self, X, y) : 
        """
        Parameters
        ----------
            y : pd.Series, these are the data to predict.  
            X : These are the features used to predict y. The goal
                of this class is to eliminate the not so useful features
                in X, to make better predictions on y. 
                transformation : transformation of y. 
        """
        self.X = X
        self.y = y        
        self.feature_names = X.columns.values 
    
    def subset_features(self, max_iter = 100000, tol = 0.001, normalize=False, inspect=False) : 
        """
        First layer of lasso regression to be used for variable selection. We keep 
        track of what features are selected to have non-zero coefs where the 
        lassoCV is fit to each year left out of the data, a bit like external 
        cross validation. 

        Parameters
        ----------
            max_iter, tol, normalize : see LassoCV() 

        """
        
        y = np.array(self.y) # needs to be np.array for maths
        X = self.X           # will be made a np.array later

        if type(X) is pd.DataFrame : 

            example_labels = X.index.values
            n_examples = X.shape[0]
            n_features = X.shape[1]
            feature_names = X.columns
            
        # Handle the type of the data 
        #if type(X) is np.ndarray :
        #    n_examples = X.shape[0]
        #    n_features = X.shape[1]
        #    example_labels = np.linspace(1, n_examples, n_examples)
        #    feature_names  = np.linspace(1, n_features, n_features)

        else :

            raise TypeError("X type: '" + str( type(X) ) + "' is not accepted. pd.Dataframe expected at this time. ")
            
        # X needs to be a np.array for the rest of this method
        X = np.array(X)
            
        if inspect :
            print(n_examples)
            print(example_labels)
            print(feature_names)
            print(n_features)
        
        # save the LassoCV alpha value for each year left out
        outer_alphas = [] 

        # Keeps track of each time a coef had a non-zero coef
        coef_counter = pd.DataFrame(0, index=[1], columns=feature_names)
        
        # for a given year left out of outside loop, what was the coef,
        # consisency will be interesting. 
        coef_values = pd.DataFrame(0, index=example_labels, columns=feature_names)
        
        # Outer loop, exclude every example once. The example iterator is what is
        # left out. 
        warnings.filterwarnings('ignore')
        count = 0
        for example in example_labels : 
            
            if inspect :
                print("leaving out: " + str(example))
                print(examples_to_keep)
            
            examples_to_keep = example_labels != example
            X_year_out       = X[examples_to_keep, :]
            y_year_out       = y[examples_to_keep]
            
            # Fit the data using leave one out 
            m = LassoCV(cv=LeaveOneOut(), normalize=normalize, tol=tol, max_iter=max_iter)
            m.fit(X_year_out, y_year_out)
            
            # TODO: Make predictions on year left out? 

            # Store fit information
            outer_alphas.append(m.alpha_)
            coef_values.loc[example, :] = m.coef_
            coef_counter.loc[1, (np.abs(m.coef_) > 0) ] += 1
            
            # progress bar
            if inspect : 
                count += 1
                print("percent completed: %f" %(100.*count/len(example_labels)))
            
        # Turn warnings back on
        warnings.filterwarnings('default')
        
        # Store these dataframes
        self.coef_counter = coef_counter
        self.coef_values  = coef_values
        self.outer_alphas = outer_alphas

        return coef_values, coef_counter, outer_alphas 

    def second_lasso(self, cutoff_percentile=0, max_iter = 100000, tol = 0.001, normalize=False) : 
        """
        Implements LassoCV() with leave-one-out to 
        find the final model. This uses the subset 
        of X values that were determined to be useful
        with the external nest
            Parameters
            ----------
                cutoff_percentile : Percentile to which to cuttoff a feature. If a feature
                                    was not selected in more than cutoff_percentile of the 
                                    external loops of subset_features() then do not include
                                    it for fitting here. When exactly zero, no feature selection
                                    is performed and all of X variables are used. 
                max_iter : see LassoCV()
                tol : see LassoCV()
                normalize: see LassoCV(), set to false because our data are centered
                           by design. 
                           
            return 
            ------
                lasso_model : The Lasso model fit using the columns of X specified as True
                              by the original_features_kept_mask. 
                X : pd.DataFrame, the dataframe used for the second fit. 
                original_features_kept_mask : Logical, the features (columbns) used (where True) 
                                     for this lasso.  
        """
        
        y = self.y
        if cutoff_percentile == 0 :
            # DO NOT subset X at all. 
            X = self.X
            print("X is not being subset")
            self.cutoff = 0
            original_features_kept_mask = np.array([True]*X.shape[1])
            
        else :
            # Subset using the specified percentile
            X_all = self.X 
             
            cutoff = np.percentile(self.coef_counter.loc[1, :].values, cutoff_percentile)
            self.cutoff = cutoff 
            original_features_kept_mask = np.array(self.coef_counter.loc[1, :].values >= cutoff)
            X = X_all.loc[:, original_features_kept_mask].copy()
            
        self.X_second_lasso = X
            
        # Store the cutoff percentile 
        self.cutoff_percentile = cutoff_percentile
        
        lasso_model = LassoCV(cv=LeaveOneOut(), normalize=normalize, tol=tol, max_iter=max_iter)
        lasso_model.fit(X.values, y)
        
        # Make this model available externally
        self.lasso_model = lasso_model 
        self.original_features_kept_mask = original_features_kept_mask
        
        # TODO: translate these coefs back to original data locaitons (columns)
        return lasso_model, X, original_features_kept_mask
    
    def plot_subset_stats(self, include_second = False, cutoff_percentile=50) :
        """
        Visualizes the results of the feature selection performed by subset_features.
        
        Parameters
        ----------
            include_second : Logical, whether to plot the results of the second lasso or not. 
            cutoff_percentile : The percentile by which features are excluded from the 
                                second use of lasso. The figures in this method visualize
                                where the percentile falls along the metadata associated 
                                with the external loop. Overidden when inclide_seconf is True. 
            
        """
        
        counts = self.coef_counter.loc[1, :].values
        
        if include_second : 
            cutoff  = self.cutoff
            cutoff_percentile = self.cutoff_percentile
        else :
            cutoff = np.percentile(counts, cutoff_percentile)
        
        n_loops = self.coef_values.shape[0]
        
        fig = plt.figure(dpi=300, figsize=(15,10))
        
        # Feature occurrence hist ----------
        ax1 = plt.subplot(231)
        h=plt.hist(counts)
        ax1.axvline(cutoff, c="r", label="cutoff for second lasso to see")
        plt.xlabel("Number of times a feature was used")
        plt.ylabel("Count for frequency")
        plt.legend()
 
        # Feature occurrence barplot ----------
        bar_colors = np.array(["C0"]*self.coef_counter.shape[1])
        if include_second :
            # Translate lasso_model non-zero coef features back
            # to full dimension of features that were used. 
            m1 = self.original_features_kept_mask 
            m2 = np.where(m1)[0][self.lasso_model.coef_!=0]
            bar_colors[m2] = "r"
        
        ax2=plt.subplot(232)
        self.coef_counter.loc[1, :].plot.bar(ax=ax2, color=bar_colors)
        ax2.axhline(cutoff, c="r")
        plt.title("feature selected out out of %i" %n_loops)
        
        # feature mean coef vs. feature occurrence scatter ----------
        ax3=plt.subplot(233)
        mean_coef_value  = np.abs(self.coef_values.mean() )
        used_count = self.coef_counter.loc[1, :] 
        plt.scatter( mean_coef_value, used_count, c = bar_colors)
        plt.xlabel("mean coef outside loop")
        plt.ylabel("count used outside loop")
        plt.title("feature mean coef vs times used")
        plt.axhline(cutoff, c="r")
        
        legend_elements = [Line2D([0], [0], marker='o', color='w', 
                                  label='original feature', markerfacecolor='C0', markersize=15),
                           Line2D([0], [0], marker='o', color='w', 
                                  label='selected by second lasso', markerfacecolor='r', markersize=15)
                          ]
        
        ax3.legend(handles=legend_elements, loc='best')
        
        # alpha distribution ----------
        ax4=plt.subplot(234)
        h2=plt.hist(self.outer_alphas)
        if include_second :
            plt.axvline(self.lasso_model.alpha_ , c="purple", label="second lasso selected alpha")
            plt.legend()
        plt.xlabel("LassoCV(leaveOneout).alpha_")
        plt.ylabel("occurrence over external loop")
        
        plt.subplots_adjust(hspace=0.7, wspace=0.3)
        
        # heat map of the coeficients from external loop ----------
        ax5 = plt.subplot(235)
        with sns.plotting_context("poster") :
            
            sns.heatmap(self.coef_values, center=0, cmap="bwr", ax=ax5)
            ax5.set(xlabel="", ylabel="year left out", title="Regression Coeficients")