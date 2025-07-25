import gzip, sys
import subprocess
import multiprocessing
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import yaml


stats = sys.argv[1].split(' ')
dynamic_filter_config = sys.argv[2]
stat_plots = sys.argv[3]
filter_combination_plot = sys.argv[4]
threads = int(sys.argv[5])
chromosomes = sys.argv[6].split(' ')
config_file = sys.argv[7]
plot_cutoff_quantile = float(sys.argv[8])
plot_cols = int(sys.argv[9])
min_filter_display_frac = float(sys.argv[10])
config_filter_sets = sys.argv[11]
filter_id = sys.argv[12]
stat_type = sys.argv[13]


with open(config_file, 'r') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
filterset = config[config_filter_sets][filter_id][stat_type]

if type(stats) == list:
    first_stats_file = stats[0]
else: first_stats_file = stats
with gzip.open(first_stats_file) as f:
    header = f.readline().decode('utf-8').strip().split()


def get_subset(fn, n_sites=None, fraction=None, header=None):
            """
            return subset dataframe
            """
            if n_sites is not None:
                # shuf: shuffle the lines; -n: number of randomly picked lines
                subset_command =  f"shuf -n {n_sites}"
            elif fraction is not None:
                subset_command =  f"awk -v seed=$RANDOM  'BEGIN{{srand(seed)}} {{if (rand()<{fraction}) print $0}}' "
            else:
                raise ValueError("n_sites or fraction must be given")
            # shell command to view stats file and get a random subset of the data
            command = f"zcat {fn} | " + subset_command
            p = subprocess.Popen([command],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, shell=True,
                                        encoding='utf8')
            s_df = pd.read_csv(p.stdout, sep='\t',na_values=['.','NA','NaN','nan'],
                                  index_col=[0,1], header=None, names=header)
            s_df = s_df.drop_duplicates(keep='first')
            o,e = p.communicate() # define errors if they occur
            if p.returncode:
                raise UserException(e)
            return s_df


def get_nsites(fn, n_sites, header=None):
    return get_subset(fn, n_sites=n_sites, fraction=None, header=header)

def get_fraction(fn, fraction, header=None):
    return get_subset(fn, n_sites=None, fraction=fraction, header=header)


if filterset['distribution_inference_subset_type'] == 'absolute':
    nsites = int(float(filterset['distribution_inference_subset']))
    frac_chrom = np.random.multinomial(nsites,(chrom_length[chromosomes]/chrom_length[chromosomes].sum()).values,)
    get_sub = get_nsites
elif filterset['distribution_inference_subset_type'] == 'fraction':
    frac = float(filterset['distribution_inference_subset'])
    frac_chrom = [frac]*len(chromosomes)
    get_sub = get_fraction
else:
    raise ValueError("distribution_inference_subset_type must be 'absolute' or 'fraction'")

    
pool = multiprocessing.Pool(threads)
procs = []
for fn, frc in zip(stats, frac_chrom):
    procs.append(pool.apply_async(get_sub, (fn, frc, header)))


stat_df = pd.concat([p.get() for p in procs])
if stat_type == 'snp':
    # keep only biallelic sites
    stat_df = stat_df[(stat_df['REF'].apply(len)==1)&(stat_df['ALT'].apply(len)==1)]
stat_df = stat_df[~stat_df.index.duplicated(keep="first")]

filters = filterset['filters'].items()


######## parse over filters ###########


def is_float(x):
    try:
        float(x)
        return True
    except:
        return False


dynamic_filter_dic = {}
filter_dic = {}


fig,axs = plt.subplots(1, plot_cols,figsize= (24,5) )
for i, (k,d) in enumerate(filters):
    tag = d['tag']
    ss = stat_df.loc[:,tag]
    ss = ss[ss.apply(is_float)].astype(float)

    try:
        d['threshold_type']
    except KeyError:
        raise KeyError("Config must provide threshold_type for filter {}.".format(k))
        
    if d['threshold_type'] == 'quantile':
        threshold = float(ss.quantile(d['threshold']))
        d.update({'threshold':threshold,
                    'threshold_type':'absolute'})
        dynamic_filter_dic.update({k:d})
    elif d['threshold_type'] == 'fraction':
        try:
            threshold = float(eval(d['summary_fun'])(ss) * d['threshold'])
            d.update({'threshold':threshold,
                    'threshold_type':'absolute'})
            dynamic_filter_dic.update({k:d})
        except KeyError:
            raise KeyError("filter {}: If threshold_type == 'fraction', "
                            "config must provide 'summary_fun'.".format(k))
    elif d['threshold_type'] == 'absolute':
        threshold = d['threshold']
    else:
        raise ValueError("Unknown threshold type {}. "
                        "Must be one of 'quantile', 'fraction', 'absolute'.".format(d['threshold_type']))        
    
    #generate dic of filter combinations
    sel_str = "ss {} {}".format(d['operator'], d['threshold'])
    filter_dic[k] = eval(sel_str)
    
    ##plot individual filters
    
    # try:
    mx = ss.quantile(plot_cutoff_quantile)
    # except TypeError:
    #     print(ss)
    #     continue
    mn = ss.min()
    ss.loc[ss > mx] = mx
    
    plt.sca(axs[i])
    ss.plot.hist(ax=axs[i], bins=100).set_xlabel(tag)
    
    if d['operator'] == '>':
        mn1 = threshold
        mx1 = mx
        pct_filtered = 100*(ss > threshold).mean()
    elif d['operator'] == '<':
        mn1 = mn
        mx1 = threshold
        pct_filtered = 100*(ss < threshold).mean()
    else:
        raise ValueError("Only < and > operators implemented for filters.")

    ymax = axs[i].get_ylim()[1]
    
    axs[i].annotate("{} {} {:.2e}".format(tag, d['operator'], threshold), (mn1 + (mx1-mn1)*0.3 , ymax*0.5))

    plt.axvspan(mn1, mx1, color='r', alpha=0.2)
    plt.title(k + " filters ~{:.2f}% of {}s".format(pct_filtered, stat_type))

plt.suptitle(filter_id + " {} filters".format(stat_type), fontsize=20)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(f"{stat_plots}.pdf")
plt.savefig(f"{stat_plots}.svg")
plt.close()


with open(dynamic_filter_config,'w') as f:
    yaml.dump(dynamic_filter_dic, f)


#evaluate filter combinations
filter_df = pd.DataFrame(filter_dic)

# total proportion filtered
filtered = filter_df.any(axis=1).mean()

filter_stats_raw = filter_df.groupby(filter_df.columns.values.tolist(), axis=0).apply(len)
filter_stats_raw.name = 'n_sites'
filter_stats = filter_stats_raw/filter_stats_raw.sum()

filter_stats_df = filter_stats.reset_index()
filter_stats_df = filter_stats_df[filter_stats_df.iloc[:,:-1].any(1)]
filter_stats_df = filter_stats_df.apply(lambda s: s.replace({True:s.name, False:''}))
filter_stats_df = filter_stats_df.set_index(filter_df.columns.values.tolist())['n_sites']

filter_stats_min_frac = filter_stats_df[filter_stats_df > min_filter_display_frac]


fig2 = plt.figure(figsize=(2*len(filter_stats_min_frac), 6))
ax = plt.gca()
(filter_stats_min_frac.sort_values(ascending=False)*100.).plot(kind='bar')
plt.title("{} filters (total {:.1f}% of {}s filtered)".format(stat_type, filtered*100, stat_type))

ax.set_ylabel("Percent filtered")
ax.set_xticklabels([", ".join([f for f in i if f]) for i in filter_stats_min_frac.index])
ax.set_xlabel("Filter combinations (>{}%)".format(100 * min_filter_display_frac))
plt.tight_layout()
plt.savefig(f"{filter_combination_plot}.pdf")
plt.savefig(f"{filter_combination_plot}.svg")
plt.close()
