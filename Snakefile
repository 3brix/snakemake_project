configfile: "config.yaml"

rule all:
    input:
        config["outputs"]["umap_plot"],
        config["outputs"]["cluster_markers"]

rule setup_and_initialization:
    input:
        data_dir=config["data_dir"]
    output:
        config["outputs"]["initial"]
    params:
        script=config["scripts"]["setup"]
    shell:
        """
        Rscript {params.script} {input.data_dir} {output}
        """

rule qc_and_preprocessing:
    input:
        config["outputs"]["initial"]
    output:
        config["outputs"]["preprocessed"]
    params:
        script=config["scripts"]["qc"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule normalization_and_feature_selection:
    input:
        config["outputs"]["preprocessed"]
    output:
        config["outputs"]["normalized"]
    params:
        script=config["scripts"]["normalization"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule scaling_and_pca:
    input:
        config["outputs"]["normalized"]
    output:
        config["outputs"]["pca"]
    params:
        script=config["scripts"]["scaling"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule clustering_and_umap:
    input:
        config["outputs"]["pca"]
    output:
        config["outputs"]["umap"]
    params:
        script=config["scripts"]["clustering"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule differential_markers_and_visualization:
    input:
        config["outputs"]["umap"]
    output:
        markers=config["outputs"]["cluster_markers"],
        plot=config["outputs"]["umap_plot_id"]
    params:
        script=config["scripts"]["markers"]
    shell:
        """
        Rscript {params.script} {input} {output.markers} {output.plot}
        """
