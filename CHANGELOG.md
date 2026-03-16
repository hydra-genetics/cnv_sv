# Changelog

## [3.0.0](https://github.com/hydra-genetics/cnv_sv/compare/cnv_sv-v2.0.0...cnv_sv-v3.0.0) (2026-03-16)


### ⚠ BREAKING CHANGES

* **sawfish:** update sawfish rules for sawfish v2 changes

### Features

* add  melt vcf post-processing ([cc4b3bb](https://github.com/hydra-genetics/cnv_sv/commit/cc4b3bb96651aad6e6b7b27d0e3faca4c2fecc50))
* add boolean vcf config ([945b44b](https://github.com/hydra-genetics/cnv_sv/commit/945b44ba8f62d2ac214e493cb2317fabcd804acf))
* add chr to chromosome if missing ([0197bb2](https://github.com/hydra-genetics/cnv_sv/commit/0197bb2c48b9d45d72c9d562a18ac02169511508))
* add entries for severus_t_only and severus_tn ([89ec119](https://github.com/hydra-genetics/cnv_sv/commit/89ec119182c830a3280a5db585abd80468ebc883))
* add function to compile input BAM paths ([49c2f1c](https://github.com/hydra-genetics/cnv_sv/commit/49c2f1c15452863db3fa9990a1cff8f7d98491ca))
* add input helper function to severus_tn; change shell commands; add bam index as input ([#261](https://github.com/hydra-genetics/cnv_sv/issues/261)) ([83b891b](https://github.com/hydra-genetics/cnv_sv/commit/83b891bf96230d7b0669211f0a0c5c97bceb8ac5))
* add melt ([6739845](https://github.com/hydra-genetics/cnv_sv/commit/6739845a9a64fc8ad0a1aab40f5e399000adbd95))
* add melt integration test files ([4710679](https://github.com/hydra-genetics/cnv_sv/commit/47106799f2e4b838c00f18b7821267197557a863))
* add melt rule and updating files connected to it ([3f3292a](https://github.com/hydra-genetics/cnv_sv/commit/3f3292a29efd38b492efd55d33621759d7c05c57))
* add melt vcf post-processing ([b11027f](https://github.com/hydra-genetics/cnv_sv/commit/b11027f6c608a53a83dd6c9170b1b793aad7b02e))
* add names for the vcfs using params and add bnd_distance param ([2f1a7f3](https://github.com/hydra-genetics/cnv_sv/commit/2f1a7f3f35f734ff5bce9e82c190c1fd2e330c57))
* add names for the vcfs using params and add bnd_distance param ([e32ca47](https://github.com/hydra-genetics/cnv_sv/commit/e32ca47bfbac69ef37bf670c9a014937c0cfd037))
* add paraphase ([6012039](https://github.com/hydra-genetics/cnv_sv/commit/6012039c4b6c90afabdc06dd82dddc7c42f09b3d))
* add paraphase ([2d533ab](https://github.com/hydra-genetics/cnv_sv/commit/2d533ab7eda83cc077803e9751f33ccae60ac05c))
* add pbsv and HiFiCNV ([c730f66](https://github.com/hydra-genetics/cnv_sv/commit/c730f66293e1ec5b8755659856347742476394af))
* add pbsv.smk ([3854216](https://github.com/hydra-genetics/cnv_sv/commit/385421622d3a7a7590b9b33fbc80515a55f3e322))
* add possibility to change caller annotation ([267285f](https://github.com/hydra-genetics/cnv_sv/commit/267285f782855a7903c4165ef829c515d2082aa2))
* add possibility to change caller name ([0b5f000](https://github.com/hydra-genetics/cnv_sv/commit/0b5f000f761b934d4ffc747b432f00efdffdebec))
* add priority as a param ([db6e1ba](https://github.com/hydra-genetics/cnv_sv/commit/db6e1bac37321f23a7d2763b94493ebc0f4dd049))
* add rule for HiFiCNV ([930b800](https://github.com/hydra-genetics/cnv_sv/commit/930b8008ea60cc36a3d9fe45d34018909aa08442))
* add savana ([f017c7d](https://github.com/hydra-genetics/cnv_sv/commit/f017c7ddce8c0723df88e33879d987fb4f0e9163))
* add savana ([1805178](https://github.com/hydra-genetics/cnv_sv/commit/1805178628156c3cd137c4df14d8396e2e355191))
* add sawfish ([6491c51](https://github.com/hydra-genetics/cnv_sv/commit/6491c51e1896adb4d3dfe116b2ef5eaf0263c493))
* add sawfish ([e7408b8](https://github.com/hydra-genetics/cnv_sv/commit/e7408b8d30c29613e389123c5fe4c87dbfa9637a))
* add scramble analysis for mobile elements ([9f5b29e](https://github.com/hydra-genetics/cnv_sv/commit/9f5b29ea109165184e3d9faf7392dd18cb5928ae))
* add scramble config ([2f92f9e](https://github.com/hydra-genetics/cnv_sv/commit/2f92f9e3d5ebcc7bc59122bdd18cd2bf44216f1a))
* add scramble in schemas ([39225c1](https://github.com/hydra-genetics/cnv_sv/commit/39225c1d9d05ebdf417859af7f7990bd8308a585))
* add scramble module ([1512273](https://github.com/hydra-genetics/cnv_sv/commit/15122738fdda141c55f55cf53c293c0ea5056cab))
* add scramble to integration config and common.smk ([2bf2f04](https://github.com/hydra-genetics/cnv_sv/commit/2bf2f047b757c24b50a52faa5870aaf703be45b4))
* add scramble to vcf rule ([49d117d](https://github.com/hydra-genetics/cnv_sv/commit/49d117d65c34b2554960d8b8a40543cdcdf8a348))
* add Severus ([7d0b7e4](https://github.com/hydra-genetics/cnv_sv/commit/7d0b7e4547f9b36e34e6c0046f8ec125adaaf360))
* add severus related entries to config.yaml ([cd64c98](https://github.com/hydra-genetics/cnv_sv/commit/cd64c9802e938117f41e2233a31a363164611f31))
* add severus.smk ([53c64dc](https://github.com/hydra-genetics/cnv_sv/commit/53c64dc0924901482bd85cb1f6f3b39dbc6b9137))
* add severus.smk to the list of included rules ([39a0dd0](https://github.com/hydra-genetics/cnv_sv/commit/39a0dd0208f1f9bf33b3b99a899880e9a091e869))
* add sniffles2 single sample calling ([1c3c8de](https://github.com/hydra-genetics/cnv_sv/commit/1c3c8de712c278ac1dbf5f2b614bfc187fbc2b9d))
* add sniffles2 single sample calling ([0ce9b8a](https://github.com/hydra-genetics/cnv_sv/commit/0ce9b8aa3821e4c98d78d8f89ca8d423164e1765))
* add tabix and bgzip for vcf output ([661268c](https://github.com/hydra-genetics/cnv_sv/commit/661268c8b8d2eafb9c6274e73e6df3b1a11932ec))
* add tabix and bgzip for vcf output ([f29e4a2](https://github.com/hydra-genetics/cnv_sv/commit/f29e4a29c469b931137f64a9698ccdb3c786c744))
* Add the CNV caller Jumble ([103c118](https://github.com/hydra-genetics/cnv_sv/commit/103c118f805830faac7ed7f93f2e9c5f7a9d5823))
* add trgt ([55a652e](https://github.com/hydra-genetics/cnv_sv/commit/55a652ece7cc9fdd4f38617963aa2703744feaf9))
* add trgt ([145641d](https://github.com/hydra-genetics/cnv_sv/commit/145641daefdd6e47d9dc07b229c964582f8874d6))
* add upd ([21c0f4d](https://github.com/hydra-genetics/cnv_sv/commit/21c0f4db7b2c48828132f662ce141e58901c0b4a))
* add upd ([cc5a660](https://github.com/hydra-genetics/cnv_sv/commit/cc5a660a61d0f6d9b4a8a4bc6559cadfcf26f741))
* added cnvkit call of jumble output ([605ce90](https://github.com/hydra-genetics/cnv_sv/commit/605ce9040f4137aead0a1f392d34cc0d931ee5fe))
* added jumble_vcf rule for segment exports ([0aa2632](https://github.com/hydra-genetics/cnv_sv/commit/0aa263211e11c6f16c56b1adb47da832cb87b905))
* added scanITD ([724e179](https://github.com/hydra-genetics/cnv_sv/commit/724e179ab2c2e958a1ce9f5a7bb097af410c5f0a))
* added scanITD ([6b7e9bf](https://github.com/hydra-genetics/cnv_sv/commit/6b7e9bf5a6ecf48c518f4542f9a6bfb2db383661))
* added the CNV caller Jumble ([4d24aa8](https://github.com/hydra-genetics/cnv_sv/commit/4d24aa883305ac8783af0df91f9b9e02e2071545))
* create empty baf, log2, probes, depth if they do not exist ([46b0939](https://github.com/hydra-genetics/cnv_sv/commit/46b0939f5a60a95bc08e0a12affd6524163e7598))
* **exomedepth:** add checkpoint rule to allow reference file as input ([a75d1e5](https://github.com/hydra-genetics/cnv_sv/commit/a75d1e5939e994ba44e665d8f6ae537228e71b5c))
* include hificnv.smk and pbsv.smk ([76ed3fc](https://github.com/hydra-genetics/cnv_sv/commit/76ed3fc59f46b5cc661b97c06f1e67300bcd2203))
* make caller name configurable ([461aee1](https://github.com/hydra-genetics/cnv_sv/commit/461aee1d4d14bc2cf337c8c634eebe2148827325))
* make the tumor cell content optional ([#160](https://github.com/hydra-genetics/cnv_sv/issues/160)) ([1c960c6](https://github.com/hydra-genetics/cnv_sv/commit/1c960c638ffdefea8b23296f8a49fbe65b82b557))
* make three mobile element rules (e.g. MELT, Scramble and xtea) ([011f900](https://github.com/hydra-genetics/cnv_sv/commit/011f900c507d874efd71b1d77989407402bb6876))
* New jumble docker ([030905f](https://github.com/hydra-genetics/cnv_sv/commit/030905f9de9011bf5b4d3f6d1c4cfb7e30184b33))
* New Jumble version with additional GIS score and changes in output names ([8007f6e](https://github.com/hydra-genetics/cnv_sv/commit/8007f6e728b99c806a547a84c05fd2af905a9dfd))
* new name and functionality of the functon for creating input bam files ([2b21dec](https://github.com/hydra-genetics/cnv_sv/commit/2b21dece2f6c8c38d3d336258c4a70793327e243))
* remove conda support and testing ([1415c8b](https://github.com/hydra-genetics/cnv_sv/commit/1415c8b5a77a68e3598ae1774e44a3d0104f67b2))
* **sawfish:** update sawfish rules for sawfish v2 changes ([2c9868f](https://github.com/hydra-genetics/cnv_sv/commit/2c9868f99cf61ce2347482913c0293cc855539a4))
* **scramble:** add unit testing ([956c1cd](https://github.com/hydra-genetics/cnv_sv/commit/956c1cda9a5d52605d50a87fe0a12bcdcece3cd4))
* update severus version to 1.5 ([e227d4f](https://github.com/hydra-genetics/cnv_sv/commit/e227d4f37155ecf6f7e3fb09da472b3b92f16659))
* use get_input_bam() in rule cnvkit_batch ([3fbd579](https://github.com/hydra-genetics/cnv_sv/commit/3fbd57986894cc5c4c917838a874e2be881d6091))
* use get_input_bam() in rule cnvkit_batch ([db571ba](https://github.com/hydra-genetics/cnv_sv/commit/db571ba3ed66f13b057e0abf575dd2afd89e8101))


### Bug Fixes

* add 'exclude' entry and 'exclude' file for hificnv ([f9c9bf5](https://github.com/hydra-genetics/cnv_sv/commit/f9c9bf50cba19fc70a08bf4ae5c39e4a7d452a45))
* add {type} to input, output, log and benchmark files in pbsv_call ([6c5c78c](https://github.com/hydra-genetics/cnv_sv/commit/6c5c78cf1b7417988114d65e8b24db7e9dad6944))
* add {type} to ouput, log and benchmark files in pbsv_discover ([f4ca49c](https://github.com/hydra-genetics/cnv_sv/commit/f4ca49c69802b0eb7656e1ffa90eab66eba417de))
* add {type} to output, log and benchmark files in rule hificnv ([26b5c38](https://github.com/hydra-genetics/cnv_sv/commit/26b5c38a4aac977af85a3b335389fabb95d75abf))
* add {unit_type} to output files of pbsv_call, pbsv_discover and hificnv in compile_output_list() ([71b2572](https://github.com/hydra-genetics/cnv_sv/commit/71b2572f902ae8db0a54cb87b5eb1b606dcb825d))
* add a config entry and a reference file required for integration test ([d1403de](https://github.com/hydra-genetics/cnv_sv/commit/d1403de4e3f0bada08fbf64edfdf6d8315ee8232))
* add all output files description to severus_* ([8705dce](https://github.com/hydra-genetics/cnv_sv/commit/8705dce34522631f47433085744f25be922e8428))
* add all required sub-entries to severus entries ([ee36340](https://github.com/hydra-genetics/cnv_sv/commit/ee36340ef96d5fcce60590c001824e86b17d82d8))
* add container for hificnv to .tests/config_pacbio.yaml ([c1fa21a](https://github.com/hydra-genetics/cnv_sv/commit/c1fa21a1b4633dc4970921fb474954a649d7ec47))
* add default values for pacbio_alignment and ont_alignment in con… ([d2a38e3](https://github.com/hydra-genetics/cnv_sv/commit/d2a38e361550a335600ef41a50759f54b96bdb50))
* add default values for pacbio_alignment and ont_alignment in config.schema.yaml ([ddedade](https://github.com/hydra-genetics/cnv_sv/commit/ddedade1728750d9cdd64f18a4d53532a05457ac))
* add dir entry to savana_* rules ([6ff1ec0](https://github.com/hydra-genetics/cnv_sv/commit/6ff1ec03dc7a9c8e93e9bcf9cadc54b8ad0fa102))
* add docker image path, reference and trf files to pbsv rules ([4842a93](https://github.com/hydra-genetics/cnv_sv/commit/4842a93e128db7f723ebea5a0138eceea81a198f))
* add entries for hificnv and pbsv ([624d28b](https://github.com/hydra-genetics/cnv_sv/commit/624d28b16292af848c290d85364c8d4b54e649d2))
* add genome version to params in the TO rules ([6c6ef84](https://github.com/hydra-genetics/cnv_sv/commit/6c6ef847a0291e5a5dbaec5e7957938ca404af0b))
* add genome version to the TO rules ([cf5da60](https://github.com/hydra-genetics/cnv_sv/commit/cf5da608e233ff744c1402e6a0a8685a18cd96ec))
* add genome_version to the TO rules ([ce3c5b6](https://github.com/hydra-genetics/cnv_sv/commit/ce3c5b6545ed2929568dc6626011a018aa672aa8))
* add imports of the new input functions to common.smk ([619b0a0](https://github.com/hydra-genetics/cnv_sv/commit/619b0a0682adea8601b3b0855134449beaf1cd5c))
* add length consensus and GT values ([a11cf27](https://github.com/hydra-genetics/cnv_sv/commit/a11cf27ec13edb3e5b1b5ed8edda9788de705daa))
* add logging ([0d5bc25](https://github.com/hydra-genetics/cnv_sv/commit/0d5bc257ad30d27107c848b023484e4463f194ad))
* add min_support to parameters ([c4f4af6](https://github.com/hydra-genetics/cnv_sv/commit/c4f4af62c77e9b343959af5eed9dc166cdd60a80))
* add missing import ([75b774c](https://github.com/hydra-genetics/cnv_sv/commit/75b774c45bd6305e656d89509c7257b6e93a085a))
* add missing out files ([9121e0c](https://github.com/hydra-genetics/cnv_sv/commit/9121e0c39f50ec18ebfc6d141dccca798091e32a))
* add missing parenthesis ([9328c59](https://github.com/hydra-genetics/cnv_sv/commit/9328c593d2bc0b3b525ef1fd8214b90e9ed2ee4a))
* add missing space in command ([21e2767](https://github.com/hydra-genetics/cnv_sv/commit/21e2767e8741a6d06cb2b4fedc5928d21025f49c))
* add output description to savana_* rules ([9dcf3d0](https://github.com/hydra-genetics/cnv_sv/commit/9dcf3d0a9225c59bb731ad0e59de36644fba9ac2))
* add output files from hificnv and pbsv to compile_output_list() ([d7241dd](https://github.com/hydra-genetics/cnv_sv/commit/d7241dd6c3a2199b07d561504e97d4f2f89bba87))
* add quite to picard for stdout ([8f09906](https://github.com/hydra-genetics/cnv_sv/commit/8f09906109194e6444e84c7b3fc6a41ecfebab73))
* add renaming of the output files after HiFiCNV run ([762fbd7](https://github.com/hydra-genetics/cnv_sv/commit/762fbd79d0f9ecc74377c585cac505e657d97d62))
* add required argument to the input function: config ([eb0eebe](https://github.com/hydra-genetics/cnv_sv/commit/eb0eebe4b0593814f90bb1ac3b39b1a19513cd09))
* add required quotation marks in awk command ([ccdc1fd](https://github.com/hydra-genetics/cnv_sv/commit/ccdc1fd373c3bb5ac3b709b0e73582d61f3f366d))
* add ruleorder directive ([9d50970](https://github.com/hydra-genetics/cnv_sv/commit/9d5097009e9a1fe627da9ff70c8fb4588e309910))
* add safer line parsing ([c921415](https://github.com/hydra-genetics/cnv_sv/commit/c921415e56c5a67791d304b651a358fbd197e8bd))
* add savana_* entries to config.yaml ([98f75e0](https://github.com/hydra-genetics/cnv_sv/commit/98f75e09332a07d83749e7d12ae98339fc0e2374))
* add scramble_vcf to schemas and softwares ([67d928c](https://github.com/hydra-genetics/cnv_sv/commit/67d928c97462ec78981bf2841d65c02b6f5e99cd))
* add temp() directive to output ([0c4add0](https://github.com/hydra-genetics/cnv_sv/commit/0c4add03f30a1322a1efdc41a29e3dc585aa9774))
* add write empty file ([4cb8450](https://github.com/hydra-genetics/cnv_sv/commit/4cb8450d0ca6917a60939e7a49437d267be50d56))
* added RTD for all reviewer rules and fix rule name ([3a81357](https://github.com/hydra-genetics/cnv_sv/commit/3a81357f9ce3d063d9a77403e06d1ea6a557b26a))
* allow missing TC when generating GATK VCF ([#162](https://github.com/hydra-genetics/cnv_sv/issues/162)) ([aa6fece](https://github.com/hydra-genetics/cnv_sv/commit/aa6fece935ea8587910b273da7cbf6b8674bfb68))
* alphabetical order of the savana_* rules in rules.schema.yaml ([238f9ce](https://github.com/hydra-genetics/cnv_sv/commit/238f9ce6fcb564aeb533e75fd3892488a4f38b25))
* automap input and integration test ([f7ed78b](https://github.com/hydra-genetics/cnv_sv/commit/f7ed78bbaba8a18292990ebb408533f3b71ca9c7))
* automap input and integration test ([20ec062](https://github.com/hydra-genetics/cnv_sv/commit/20ec062be15dce5a3766f50254d3fdd73b86361e))
* better handle sample.tsv without a sex column ([90696c5](https://github.com/hydra-genetics/cnv_sv/commit/90696c5d35715376755c8ee35a978bdc55fc0b56))
* better handle sample.tsv without a sex column ([990acf6](https://github.com/hydra-genetics/cnv_sv/commit/990acf698ce6a8613b1217c4dd23ee31d49bac0c))
* bump snakemake version to support apptainer ([#151](https://github.com/hydra-genetics/cnv_sv/issues/151)) ([eb4b05d](https://github.com/hydra-genetics/cnv_sv/commit/eb4b05df754c2b6be16d5feed3244d00d9f4b5f1))
* change input BAM files in pbsv_discover from hard-coded to function ([cab49eb](https://github.com/hydra-genetics/cnv_sv/commit/cab49ebb15aecac2e6b53d62b4af7b333b54f00c))
* change input bam path; mv ref to params ([6075a12](https://github.com/hydra-genetics/cnv_sv/commit/6075a12fdfd8d355c1ee391b8168616ec9d9fe7f))
* change input from  path to function ([8b7c139](https://github.com/hydra-genetics/cnv_sv/commit/8b7c139bff70774c2953654ba559c72d7a16bb6f))
* change output list gathering ([80983c2](https://github.com/hydra-genetics/cnv_sv/commit/80983c2e3c5a203f064d36e6a24b92f7e763ed9e))
* change output/log to vcf.gz ([7d409db](https://github.com/hydra-genetics/cnv_sv/commit/7d409dbf66e7e74cf6204e8e3a14afa9623dd47a))
* change to compatible wrapper ([4a83905](https://github.com/hydra-genetics/cnv_sv/commit/4a83905be016d3b6bc781b89faa4bf65588f8a89))
* change to compatible wrapper ([4759cff](https://github.com/hydra-genetics/cnv_sv/commit/4759cffbb9cbd2ac24efb85296cb9a854fcc19ef))
* change to config.get in exomedepth_call ([feeb98e](https://github.com/hydra-genetics/cnv_sv/commit/feeb98e48ab922ca9d59d0badd787a34107bbd8f))
* change to config.get in exomedepth_call ([34f5753](https://github.com/hydra-genetics/cnv_sv/commit/34f5753530450e6ae8ba51bd19c5bb0946a45d14))
* cnvkit batch v9.2.0 ([#259](https://github.com/hydra-genetics/cnv_sv/issues/259)) ([d3e4b65](https://github.com/hydra-genetics/cnv_sv/commit/d3e4b65d2f088e0975eeb85f3ee89d2e3ae763b5))
* **cnvkit:** update docker version ([9d90ccb](https://github.com/hydra-genetics/cnv_sv/commit/9d90ccb4b9161980cbf0f32c2662712e142779d6))
* **cnvpytor:** add bam index as input ([f60afd3](https://github.com/hydra-genetics/cnv_sv/commit/f60afd3386076dc2eee6806c4ff914aa2a64b763))
* **cnvpytor:** add bam index as input ([60c840a](https://github.com/hydra-genetics/cnv_sv/commit/60c840ad7c71dbfd66c0925becc210ae31dfa41a))
* code improvement ([911422a](https://github.com/hydra-genetics/cnv_sv/commit/911422a2cd2428c574b7f398c8129f0a340577a6))
* **common:** set type to str for pandas dataframe ([2d76ea3](https://github.com/hydra-genetics/cnv_sv/commit/2d76ea3370e3fa324498423a521686806a6119ef))
* **common:** set type to str for pandas dataframe ([4c0e0f2](https://github.com/hydra-genetics/cnv_sv/commit/4c0e0f2252da1ab187ff94c132523fe1958c959c))
* correct alignment_path and index_path in get_input_bam() ([d232314](https://github.com/hydra-genetics/cnv_sv/commit/d2323143815b5bc19854817f8bbd4fdb4a539b38))
* correct config entry name for pon ([c867311](https://github.com/hydra-genetics/cnv_sv/commit/c867311cd4e644b78ea21c10386bb5010ea3f25b))
* correct input parameter name ([1722db4](https://github.com/hydra-genetics/cnv_sv/commit/1722db4fea07392d23eb31c13bf6bd4586fb5ae6))
* correct output subdir name for severus_tn ([74c0648](https://github.com/hydra-genetics/cnv_sv/commit/74c064801b0ab0efcb7c1a7b10e476852a36aeb0))
* corrected output file names ([238fb36](https://github.com/hydra-genetics/cnv_sv/commit/238fb36ac6e9ffe4bb6318ef244e564e5f7719ce))
* Delete workflow/rules/juli.smk ([d72fc1e](https://github.com/hydra-genetics/cnv_sv/commit/d72fc1e5dd68333144439d86cb24fe8c4abc1abe))
* **exomedepth:** chosse female ref when the peddy sex is NA ([eba758d](https://github.com/hydra-genetics/cnv_sv/commit/eba758d79e31b634142a9d22ad9ff50a83b5fa28))
* **exomedepth:** chosse female ref when the peddy sex is NA ([4357702](https://github.com/hydra-genetics/cnv_sv/commit/4357702e36f4b472b57bf6643ea19f26bf6356b0))
* full output filenames & outdir to params ([5d23ac2](https://github.com/hydra-genetics/cnv_sv/commit/5d23ac28507a1d1a7f7a8e1dba1eeb2518ec397c))
* get mobile_elements up to date ([899ad2d](https://github.com/hydra-genetics/cnv_sv/commit/899ad2d01939db24bf9e3a3a8a8b7a28799c0f29))
* handle that cns files can have different number of columns ([203e60a](https://github.com/hydra-genetics/cnv_sv/commit/203e60a9d0c49e6682b7b6a1a7e5af4c046e4559))
* handle that cns files can have different number of columns ([57da77b](https://github.com/hydra-genetics/cnv_sv/commit/57da77b993ea25be88116e75712b09411e45677e))
* **hificnv:** add type in output files from hificnv, to reflect the new addotion of type in BAM RG SM ([b7b69c4](https://github.com/hydra-genetics/cnv_sv/commit/b7b69c4a2489eb7c9db870420138720735cd8354))
* incorrect type for copy number thresholds ([#164](https://github.com/hydra-genetics/cnv_sv/issues/164)) ([b549266](https://github.com/hydra-genetics/cnv_sv/commit/b549266994149bf7b4bdcc385456229b4f0eee77))
* jumble results removed prematurely ([ec200f8](https://github.com/hydra-genetics/cnv_sv/commit/ec200f89c9bd0f661ac391f23f05175cda6a828f))
* jumble results removed prematurely ([5ae3b74](https://github.com/hydra-genetics/cnv_sv/commit/5ae3b74b4da103c03332c77f2caa19aed5e2bee1))
* jumble_vcf should have same parameters as cnvkit_vcf ([1748415](https://github.com/hydra-genetics/cnv_sv/commit/1748415286965c99ea7864c5fc11694f3f8f6568))
* limit number of threads used by Jumble to that specified in config ([5f8d984](https://github.com/hydra-genetics/cnv_sv/commit/5f8d98418fe7cbb6b286666a3e0c2dcbc86a6abe))
* limit number of threads used by Jumble to that specified in config ([b27fb17](https://github.com/hydra-genetics/cnv_sv/commit/b27fb174baa146d5492e787232f6cb4a3805477f))
* linting ([e28d277](https://github.com/hydra-genetics/cnv_sv/commit/e28d27753c0373190050474ecee1b5637151f596))
* make name configurable ([69ea10a](https://github.com/hydra-genetics/cnv_sv/commit/69ea10ac806fd0fec770a72e3974d7e1b6a18a14))
* make sure sample column is a string ([2d3aa6b](https://github.com/hydra-genetics/cnv_sv/commit/2d3aa6bae5addca7219e9aa833778110efd1173f))
* make two entries for pbsv (call & discover) ([c8a9272](https://github.com/hydra-genetics/cnv_sv/commit/c8a9272e129f1bd49c9b355ada28353dc075995f))
* make two entries for pbsv (call & discover) ([55b8793](https://github.com/hydra-genetics/cnv_sv/commit/55b8793b11a3b0688a42beea99474c458251bdfa))
* **melt:** move log and benchmark out of temped directory ([5ee6370](https://github.com/hydra-genetics/cnv_sv/commit/5ee63704f741ca300a90d081958c65d1751581c7))
* **melt:** move log and benchmark out of temped directory ([9fe798a](https://github.com/hydra-genetics/cnv_sv/commit/9fe798ae44c3b4878d2348f8a09f70a2df3bc3a9))
* **melt:** remove hardcoded path ([9f77770](https://github.com/hydra-genetics/cnv_sv/commit/9f77770dd5980b39d1fa9df6ec4ea3abcd540723))
* **melt:** remove hardcoded path ([b15f928](https://github.com/hydra-genetics/cnv_sv/commit/b15f9289ae1c6fbaa88bf40bb50f864c7dfd73bf))
* **melt:** remove hardcoded path ([f441cb7](https://github.com/hydra-genetics/cnv_sv/commit/f441cb77e34d4223f96ef2d0a3e23527cf031176))
* **melt:** Update output file paths in melt.smk ([a812ae7](https://github.com/hydra-genetics/cnv_sv/commit/a812ae7a04aad2acd2c3d19e87cff51697626c2b))
* **melt:** Update output file paths in melt.smk ([37070c6](https://github.com/hydra-genetics/cnv_sv/commit/37070c6d954ff1a789d32d2e1240fc140b110feb))
* missing parenthesis ([be912f8](https://github.com/hydra-genetics/cnv_sv/commit/be912f822c7f7284676c7b083dc51214da5bb9d8))
* new input function in cnvkit batch ([9db74c1](https://github.com/hydra-genetics/cnv_sv/commit/9db74c13c95cd48ea0f715690bd8f6da54b3ce1e))
* optional bed file ([7e2d882](https://github.com/hydra-genetics/cnv_sv/commit/7e2d8821e7206859d58532f0889455f2abf01ef3))
* **paraphase:** add missing extra to shell command ([1250f24](https://github.com/hydra-genetics/cnv_sv/commit/1250f243de2623e1c235a26258ac32febceafa45))
* pin pulp version to &lt;2.8.0 ([1c960c6](https://github.com/hydra-genetics/cnv_sv/commit/1c960c638ffdefea8b23296f8a49fbe65b82b557))
* pindel vcf add correct samplename and contigs in header ([34e6a66](https://github.com/hydra-genetics/cnv_sv/commit/34e6a66f0a548a3021d4dd3f6792a72c856daff3))
* put param in correct rule ([b99712e](https://github.com/hydra-genetics/cnv_sv/commit/b99712e44513dc7840a1da9fac5b0244ecb60b15))
* remove duplicate {sample} wildcard from severus output paths ([#263](https://github.com/hydra-genetics/cnv_sv/issues/263)) ([21d3697](https://github.com/hydra-genetics/cnv_sv/commit/21d3697a459a49b171bd81c84ce0bc16f3a90fbe))
* remove import get_longread_bam ([cb26475](https://github.com/hydra-genetics/cnv_sv/commit/cb26475c2062d46441ec7b215ba9e7d807f8350f))
* remove import get_longread_bam ([d288e91](https://github.com/hydra-genetics/cnv_sv/commit/d288e91d08458afe7714a109b70c2c99675af8f9))
* remove print() from get_input_bam() ([80800a4](https://github.com/hydra-genetics/cnv_sv/commit/80800a4e1724d0920f6a6ea7a1487f2ec3fa57f7))
* remove tmp from final output ([0179e8b](https://github.com/hydra-genetics/cnv_sv/commit/0179e8b7c05430fb556c8dc7adfd6323079ab55c))
* replace the input funcitno in rule cnvkit_batch ([da27dc4](https://github.com/hydra-genetics/cnv_sv/commit/da27dc4e61333629083b2ad158171db50faeb865))
* rm extra } in benchmark name ([b4d01a9](https://github.com/hydra-genetics/cnv_sv/commit/b4d01a93398b1d96d58035138567c2b6545b7720))
* rm extra } in benchmark name ([e2f34ce](https://github.com/hydra-genetics/cnv_sv/commit/e2f34ce8d24515d7d1f980994f3d72e12878d4da))
* **sawfish:** update inputs to function ([d8742a6](https://github.com/hydra-genetics/cnv_sv/commit/d8742a647e50e4b6f70fda530713a6199750b6db))
* **sawfish:** update inputs to function ([671b316](https://github.com/hydra-genetics/cnv_sv/commit/671b316e8f2c7410c08447f0ae959e3862bbd9c0))
* **scanITD:** correct rule params ([01cb157](https://github.com/hydra-genetics/cnv_sv/commit/01cb157b762dd03023eea58aefcadd74a2fde9a6))
* **scanITD:** correct shell command ([3da1a76](https://github.com/hydra-genetics/cnv_sv/commit/3da1a76bf95766149f9d750f9472d8ec3ed1b185))
* **scramble:** update config ([9ad3114](https://github.com/hydra-genetics/cnv_sv/commit/9ad311447f00363609e7fe3634364c3552ca3880))
* **scramble:** update config to match rule ([f1fbf11](https://github.com/hydra-genetics/cnv_sv/commit/f1fbf11c8a3f36f9a7bf0531d17ecd6733b0cf02))
* **scramble:** update schemas to match rule ([d00f6ab](https://github.com/hydra-genetics/cnv_sv/commit/d00f6abcd2119e88cfd734d04abc76eb6ff7d3ac))
* separate output files not only dir ([f977330](https://github.com/hydra-genetics/cnv_sv/commit/f97733055ef615689608e95680864fa528b69058))
* set sample_name correctly in header ([aa82c08](https://github.com/hydra-genetics/cnv_sv/commit/aa82c085ea8115cd6700e1e91feaa61195e42930))
* set sample_name correctly in header ([ef428c6](https://github.com/hydra-genetics/cnv_sv/commit/ef428c6d69394fd8495eaf8ebfe5db730932f475))
* skip making hard links before copy ([4d8a1f5](https://github.com/hydra-genetics/cnv_sv/commit/4d8a1f54f6e1b5f7b64281cc487320842e6d54f2))
* skip making hard links before copy ([2989636](https://github.com/hydra-genetics/cnv_sv/commit/2989636f36e92723ebbff139a88bd5f074c0f13e))
* **sniffles2:** add sample-id param to add sample name in VCF ([56b818f](https://github.com/hydra-genetics/cnv_sv/commit/56b818fcbff75ab18036c8e6a4a21c92ab8a2bc0))
* **sniffles2:** add sample-id param to add sample name in VCF ([7f236b9](https://github.com/hydra-genetics/cnv_sv/commit/7f236b990664d66de412663bc4cf4913836187bf))
* solve conditional output for vcf ([a705f97](https://github.com/hydra-genetics/cnv_sv/commit/a705f97040d55cf1e956a00a0de9dd1d976eb1ea))
* solve test data mismatch ([1a1af64](https://github.com/hydra-genetics/cnv_sv/commit/1a1af64becaefe5eb7beaa606b749af96617b314))
* solve vcf output and reference dependency ([fa8b556](https://github.com/hydra-genetics/cnv_sv/commit/fa8b5568ee9264811482c5c685bb09110c3d4ef8))
* specify all output files and fix the 'ouput file name as prefix' problem ([8d9490e](https://github.com/hydra-genetics/cnv_sv/commit/8d9490ef4c232a1c0dcf8124b22f0af79a612ae0))
* specify model when calling and filtering in cnvpytor, and thread… ([0b0cd74](https://github.com/hydra-genetics/cnv_sv/commit/0b0cd74eaa2501bcad3affc018baec63d50dd89d))
* specify model when calling and filtering in cnvpytor, and threads in readdpth rule ([4347806](https://github.com/hydra-genetics/cnv_sv/commit/434780660b693fc7b9dbdf709cc4700be33de63c))
* **svdb:** use the overlap and extra params ([7e77719](https://github.com/hydra-genetics/cnv_sv/commit/7e77719ff4f2ba16457e734154e54cc32ad5f337))
* **svdb:** use the overlap and extra params ([fe656a9](https://github.com/hydra-genetics/cnv_sv/commit/fe656a91dd0a74c566ef0b8bf3512609a071ab94))
* **tiddit:** add bai file as an input ([cd88fdb](https://github.com/hydra-genetics/cnv_sv/commit/cd88fdbc61b29149610685005f6ab055340ed971))
* **tiddit:** add bai file as an input ([ace4bcd](https://github.com/hydra-genetics/cnv_sv/commit/ace4bcd17b2f7c16f58b9eac6b223308067a9d45))
* **tiddit:** add threads and access extras ([4e0fc0c](https://github.com/hydra-genetics/cnv_sv/commit/4e0fc0c2f688e5dc234560e53f985fe7150e5ac6))
* **tiddit:** add threads and access extras ([330c506](https://github.com/hydra-genetics/cnv_sv/commit/330c5064b76edbb8205b6d630f311a3dae3ff4ca))
* update common container ([05d05dd](https://github.com/hydra-genetics/cnv_sv/commit/05d05dda89ba8b1a3eb2b5777d8d54f664c563f0))
* update common container ([409dbda](https://github.com/hydra-genetics/cnv_sv/commit/409dbda6607a2d6cbba60ea63632bebb9893e333))
* Update config_pacbio.yaml ([e7f6117](https://github.com/hydra-genetics/cnv_sv/commit/e7f6117eac6d519f1238395b0bb79864b6d9daa9))
* Update config.yaml ([c2a82f8](https://github.com/hydra-genetics/cnv_sv/commit/c2a82f81e283dabec54a7e8494ba1adf41a2b574))
* Update config.yaml ([ea9a0f7](https://github.com/hydra-genetics/cnv_sv/commit/ea9a0f7d916952197381fcb61c98063a2f71835b))
* Update config.yaml ([794e3d6](https://github.com/hydra-genetics/cnv_sv/commit/794e3d6fedc8d26674128bb7f06db814e7900d44))
* update container for rule pindel_update_vcf ([#161](https://github.com/hydra-genetics/cnv_sv/issues/161)) ([512a924](https://github.com/hydra-genetics/cnv_sv/commit/512a924beba54282757745c8c66475e78f52cca3))
* update container version ([58feca4](https://github.com/hydra-genetics/cnv_sv/commit/58feca4a8f9c728e175646f9b66a2c9880b5d6fb))
* update input / output dependencies ([6338879](https://github.com/hydra-genetics/cnv_sv/commit/633887965b04c3918cb0169e865552f714a25639))
* update to bugfixed docker ([20477be](https://github.com/hydra-genetics/cnv_sv/commit/20477be5d70c92803b258398089e478015289f6b))
* update to hydra-genetics 3.3.0 ([7887b05](https://github.com/hydra-genetics/cnv_sv/commit/7887b0575296bf5971126e32c02c1e8bf3353e53))
* use correct name of the input bam function ([0ff7ccd](https://github.com/hydra-genetics/cnv_sv/commit/0ff7ccdf1c4eb70f3990fb25bf960c848e2c7de6))
* use get_input_bam() ([5aeb9af](https://github.com/hydra-genetics/cnv_sv/commit/5aeb9af6deb282743a7191430f78bd8b63ccf85a))
* use new input funciton for aligned bam files ([c4d1b2b](https://github.com/hydra-genetics/cnv_sv/commit/c4d1b2bc677b638f8baf7f2d7fa12041ce07221c))
* wrap in main ([73a5aa7](https://github.com/hydra-genetics/cnv_sv/commit/73a5aa79d29e2318000501a7cada1387dcefad02))
* wrong config parameter name ([3ad5336](https://github.com/hydra-genetics/cnv_sv/commit/3ad533667c4ea0577793c7e55cdb3cf3474cb23e))


### Performance Improvements

* add bcftools reheader to pindel vcf for correct samplename ([fcadaba](https://github.com/hydra-genetics/cnv_sv/commit/fcadaba56976bad8a6f2975558c5368c25016fb8))


### Documentation

* add a link to the input functions page on HG read-the-docs ([d2054a3](https://github.com/hydra-genetics/cnv_sv/commit/d2054a3b63895fdfbf4df069ce438fff70c66a87))
* add description of savana_* rules to softwares.md ([a2d88c7](https://github.com/hydra-genetics/cnv_sv/commit/a2d88c790e59be2773ad173b0bf32ea0af7a28c5))
* add entries for savana_* rules to config.schema.yaml ([8a14af8](https://github.com/hydra-genetics/cnv_sv/commit/8a14af8fcc9b80837b12131e735ce89127c003de))
* add schemas for rules savana_* ([c9b1a43](https://github.com/hydra-genetics/cnv_sv/commit/c9b1a43b9df5ca13a3b542f08d4980605424039e))
* added docs and tests ([0fc6464](https://github.com/hydra-genetics/cnv_sv/commit/0fc64646e1049d93a85021ce1a890220d37582f9))
* Added read the docs for all rules! ([c55a89e](https://github.com/hydra-genetics/cnv_sv/commit/c55a89e10fc5dadf3a780d463967435903bbc238))
* added RTD files and the rules for automap and cnvkit ([ce2f2ca](https://github.com/hydra-genetics/cnv_sv/commit/ce2f2ca0b4cdd6876968f638d3d30226f62bfe19))
* added RTD for all manta rules ([2571744](https://github.com/hydra-genetics/cnv_sv/commit/25717444f3a6a48825e95d2c4e756ec2af62c244))
* added RTD for all pindel rules ([b0ca4da](https://github.com/hydra-genetics/cnv_sv/commit/b0ca4dadbbada7bbad66c1274c125f58ef69447d))
* added RTD for all purecn rules ([7211b29](https://github.com/hydra-genetics/cnv_sv/commit/7211b298c138c8aaa693bca072b78885379d8532))
* added RTD for all smncaller rules ([ed55b6b](https://github.com/hydra-genetics/cnv_sv/commit/ed55b6b121c08fb4b4935fe231f493e563f41104))
* added RTD for all svdb rules ([1cdd6ba](https://github.com/hydra-genetics/cnv_sv/commit/1cdd6ba37ca9cfd3c21458188d3d04b706334112))
* added RTD for cnvpytor and exomedepth ([8fd490a](https://github.com/hydra-genetics/cnv_sv/commit/8fd490a8b73efc412797859678569ecba07ac783))
* added RTD for expansionhunter ([4f290ca](https://github.com/hydra-genetics/cnv_sv/commit/4f290ca1e0b69dc7540411d877d380644da69b7e))
* added RTD for gatk cnv rules ([a0a575e](https://github.com/hydra-genetics/cnv_sv/commit/a0a575ef4a37ac0014b67da8690f8241eeef68f9))
* added RTD for tiddit rules ([fea5964](https://github.com/hydra-genetics/cnv_sv/commit/fea5964ec64f5e184c3df75a43241485fbdfa8b0))
* added RTD for upd rules ([2eccfe5](https://github.com/hydra-genetics/cnv_sv/commit/2eccfe5bc01c485afc3996cd85cf0c3700a06175))
* dag graph ([7a654c2](https://github.com/hydra-genetics/cnv_sv/commit/7a654c2ddedecc9a47321f8120b955162fe66674))
* documentation fix ([958550b](https://github.com/hydra-genetics/cnv_sv/commit/958550be4c0b91ac499f73310ce4f4078e76b6b0))
* fix of format string ([eec358b](https://github.com/hydra-genetics/cnv_sv/commit/eec358b09903abd3038f3b5386565d51e6542e8c))
* organise rules based on software ([36568d9](https://github.com/hydra-genetics/cnv_sv/commit/36568d9cc61789b3ee5e07af1d695c01d783e496))
* remove dp_bw as output for sawfish discover ([b661ad0](https://github.com/hydra-genetics/cnv_sv/commit/b661ad0a8dcb89ef79e1ecdf6e84db9f5061a657))
* remove rules not used ([08e3600](https://github.com/hydra-genetics/cnv_sv/commit/08e360041749bfde0cb12664ba0380adcc801e4c))
* schema updates ([e07415f](https://github.com/hydra-genetics/cnv_sv/commit/e07415fa2aa238f17fecea6ceea8b29b5072bb76))
* update compatibility ([e6011a5](https://github.com/hydra-genetics/cnv_sv/commit/e6011a539ce2d6e5d486fc2feaedbce2e1c517a1))
* update compatibility ([0e1e381](https://github.com/hydra-genetics/cnv_sv/commit/0e1e381ac649bb48a0573f932053d48153e86beb))
* update compatibility ([fe31a4b](https://github.com/hydra-genetics/cnv_sv/commit/fe31a4b4b39cb5413c8875947d22156a416ea1e2))
* update DAGs and README.md ([1b0bbd0](https://github.com/hydra-genetics/cnv_sv/commit/1b0bbd0589ba72a457491068f86ff50c135f5ed1))
* update DAGs and README.md ([47afb28](https://github.com/hydra-genetics/cnv_sv/commit/47afb28ec3545f1ab5d1f4e0460d440986c6f602))
* update disclaimer ([77e41f2](https://github.com/hydra-genetics/cnv_sv/commit/77e41f2d908eb9e8e821f160283c87c1bce86831))
* update disclaimer ([1d98fea](https://github.com/hydra-genetics/cnv_sv/commit/1d98feab3ab2f1c6be40909119c00baee9045a48))
* update disclaimer ([cbb7592](https://github.com/hydra-genetics/cnv_sv/commit/cbb7592ecf98ffb4245e7673bcbefde2a4759a82))
* update docs/softwares.md ([d5fdce8](https://github.com/hydra-genetics/cnv_sv/commit/d5fdce851f04adb90f60f8588b387e99be8024bb))
* update docs/softwares.md ([688c081](https://github.com/hydra-genetics/cnv_sv/commit/688c08145e10f4cbed7c608cf346dc52353ec07b))
* update docs/softwares.md ([d76e841](https://github.com/hydra-genetics/cnv_sv/commit/d76e8411c74cec9f73e350bd48d9eb3254050720))
* update docs/softwares.md ([4268631](https://github.com/hydra-genetics/cnv_sv/commit/4268631f82c8fef281c584ccc25a9949f22f4ade))
* Update docs/softwares.md ([d3e3f67](https://github.com/hydra-genetics/cnv_sv/commit/d3e3f67f1e91e15223cc655133c944128c9e5a54))
* Update docs/softwares.md ([188ef59](https://github.com/hydra-genetics/cnv_sv/commit/188ef594302900cc0b37580e616775412fe79764))
* Update docs/softwares.md ([90ad983](https://github.com/hydra-genetics/cnv_sv/commit/90ad983f53a88117f0532bb8042c74c279d1d7d5))
* Update docs/softwares.md ([dd02db5](https://github.com/hydra-genetics/cnv_sv/commit/dd02db57a7d41c253c1e1d5e79b521d43a6f003d))
* Update docs/softwares.md ([6f16c13](https://github.com/hydra-genetics/cnv_sv/commit/6f16c13479c785cfa52f7c318a9c85aa6a2f5727))
* Update docs/softwares.md ([25b3195](https://github.com/hydra-genetics/cnv_sv/commit/25b31959177b1cd798e4896e47c5bab36c94bad3))
* Update docs/softwares.md ([c02990c](https://github.com/hydra-genetics/cnv_sv/commit/c02990c978a76a7d26ed8f1ce4f962e42bc566e3))
* update documentation ([ddd4f8c](https://github.com/hydra-genetics/cnv_sv/commit/ddd4f8c392e13900aed2b258a9d430344839d3a2))
* update documentation ([087f67e](https://github.com/hydra-genetics/cnv_sv/commit/087f67e52ce4cda10861f58f4c42609b354dad95))
* update message format ([3ae7e19](https://github.com/hydra-genetics/cnv_sv/commit/3ae7e19b179e8ce9ef6a989e6dc8c61e33af342b))
* update message format ([76eaa0b](https://github.com/hydra-genetics/cnv_sv/commit/76eaa0bca03ee6942ee83b63dd68f5bf256f364d))
* update message format ([4e4f73a](https://github.com/hydra-genetics/cnv_sv/commit/4e4f73ad6f3c3e2ba5cd0896475ee63e8ebe5533))
* update mkdocs ([e28d589](https://github.com/hydra-genetics/cnv_sv/commit/e28d5896c9dbd3727bcff704bf6e7967b2b43bf0))
* update pacbio dag and restructure module output files ([c02dc01](https://github.com/hydra-genetics/cnv_sv/commit/c02dc01d4910c3d71687adbbe1d1982b399aa902))
* update pacbio dag and restructure module output files ([876faa8](https://github.com/hydra-genetics/cnv_sv/commit/876faa8d06f41912239bf97900ea83e99fa4e5d2))
* update README.md ([eb7afab](https://github.com/hydra-genetics/cnv_sv/commit/eb7afab1fa264a1454611bb504f849e9bcacc98a))
* update README.md ([35c026d](https://github.com/hydra-genetics/cnv_sv/commit/35c026dd0d50465bac469e725aa15fe0b1a77b5c))
* update rule name for docs ([cbd7ee2](https://github.com/hydra-genetics/cnv_sv/commit/cbd7ee218eb1c0983a7f959600bbe72a354e129b))
* update rules ([e1166af](https://github.com/hydra-genetics/cnv_sv/commit/e1166aff1bee8ace4b6f4a8583cfd9d7a4cd7218))
* update schemas with sawfish info ([e8f0163](https://github.com/hydra-genetics/cnv_sv/commit/e8f0163b632499f4f81908dc5d76d7d24c48e763))
* update the module output files ([02aca91](https://github.com/hydra-genetics/cnv_sv/commit/02aca91417bfabed6eb6ee1dc8583c609b9cc02b))
* update the rule schema and softwares.md with changes to sawfish rules ([9635bec](https://github.com/hydra-genetics/cnv_sv/commit/9635bec0749a3419f4c3a89546dcc2955ca2c3ae))
* update to new rule plugin ([aedbdfd](https://github.com/hydra-genetics/cnv_sv/commit/aedbdfd8d9f3084b4e49a0d7e9bb709da707fb10))
* Update workflow/schemas/resources.schema.yaml ([6136f71](https://github.com/hydra-genetics/cnv_sv/commit/6136f71519cb149297621d16cb16c31bf8017838))
* updated output files ([c5d3d62](https://github.com/hydra-genetics/cnv_sv/commit/c5d3d625df5a560dddb32d2d74054558e6ff7565))
* updated output files ([d56c8d6](https://github.com/hydra-genetics/cnv_sv/commit/d56c8d6282382719dc2ce0c5751a2097f616f8b6))

### [1.0.1](https://www.github.com/hydra-genetics/cnv_sv/compare/v1.0.0...v1.0.1) (2025-09-01)


### Bug Fixes

* better handle sample.tsv without a sex column ([990acf6](https://www.github.com/hydra-genetics/cnv_sv/commit/990acf698ce6a8613b1217c4dd23ee31d49bac0c))


### Documentation

* update README.md ([35c026d](https://www.github.com/hydra-genetics/cnv_sv/commit/35c026dd0d50465bac469e725aa15fe0b1a77b5c))

## [1.0.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.10.0...v1.0.0) (2025-08-29)


### ⚠ BREAKING CHANGES

* **sawfish:** update sawfish rules for sawfish v2 changes

### Features

* **sawfish:** update sawfish rules for sawfish v2 changes ([2c9868f](https://www.github.com/hydra-genetics/cnv_sv/commit/2c9868f99cf61ce2347482913c0293cc855539a4))


### Bug Fixes

* **hificnv:** add type in output files from hificnv, to reflect the new addotion of type in BAM RG SM ([b7b69c4](https://www.github.com/hydra-genetics/cnv_sv/commit/b7b69c4a2489eb7c9db870420138720735cd8354))


### Documentation

* remove dp_bw as output for sawfish discover ([b661ad0](https://www.github.com/hydra-genetics/cnv_sv/commit/b661ad0a8dcb89ef79e1ecdf6e84db9f5061a657))
* update pacbio dag and restructure module output files ([876faa8](https://www.github.com/hydra-genetics/cnv_sv/commit/876faa8d06f41912239bf97900ea83e99fa4e5d2))
* update rule name for docs ([cbd7ee2](https://www.github.com/hydra-genetics/cnv_sv/commit/cbd7ee218eb1c0983a7f959600bbe72a354e129b))
* update the rule schema and softwares.md with changes to sawfish rules ([9635bec](https://www.github.com/hydra-genetics/cnv_sv/commit/9635bec0749a3419f4c3a89546dcc2955ca2c3ae))

## [0.10.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.9.0...v0.10.0) (2025-08-27)


### Features

* add savana ([1805178](https://www.github.com/hydra-genetics/cnv_sv/commit/1805178628156c3cd137c4df14d8396e2e355191))


### Bug Fixes

* add dir entry to savana_* rules ([6ff1ec0](https://www.github.com/hydra-genetics/cnv_sv/commit/6ff1ec03dc7a9c8e93e9bcf9cadc54b8ad0fa102))
* add genome version to params in the TO rules ([6c6ef84](https://www.github.com/hydra-genetics/cnv_sv/commit/6c6ef847a0291e5a5dbaec5e7957938ca404af0b))
* add genome version to the TO rules ([cf5da60](https://www.github.com/hydra-genetics/cnv_sv/commit/cf5da608e233ff744c1402e6a0a8685a18cd96ec))
* add genome_version to the TO rules ([ce3c5b6](https://www.github.com/hydra-genetics/cnv_sv/commit/ce3c5b6545ed2929568dc6626011a018aa672aa8))
* add min_support to parameters ([c4f4af6](https://www.github.com/hydra-genetics/cnv_sv/commit/c4f4af62c77e9b343959af5eed9dc166cdd60a80))
* add output description to savana_* rules ([9dcf3d0](https://www.github.com/hydra-genetics/cnv_sv/commit/9dcf3d0a9225c59bb731ad0e59de36644fba9ac2))
* add savana_* entries to config.yaml ([98f75e0](https://www.github.com/hydra-genetics/cnv_sv/commit/98f75e09332a07d83749e7d12ae98339fc0e2374))
* alphabetical order of the savana_* rules in rules.schema.yaml ([238f9ce](https://www.github.com/hydra-genetics/cnv_sv/commit/238f9ce6fcb564aeb533e75fd3892488a4f38b25))
* full output filenames & outdir to params ([5d23ac2](https://www.github.com/hydra-genetics/cnv_sv/commit/5d23ac28507a1d1a7f7a8e1dba1eeb2518ec397c))


### Documentation

* add description of savana_* rules to softwares.md ([a2d88c7](https://www.github.com/hydra-genetics/cnv_sv/commit/a2d88c790e59be2773ad173b0bf32ea0af7a28c5))
* add entries for savana_* rules to config.schema.yaml ([8a14af8](https://www.github.com/hydra-genetics/cnv_sv/commit/8a14af8fcc9b80837b12131e735ce89127c003de))
* add schemas for rules savana_* ([c9b1a43](https://www.github.com/hydra-genetics/cnv_sv/commit/c9b1a43b9df5ca13a3b542f08d4980605424039e))
* update DAGs and README.md ([47afb28](https://www.github.com/hydra-genetics/cnv_sv/commit/47afb28ec3545f1ab5d1f4e0460d440986c6f602))
* update docs/softwares.md ([d5fdce8](https://www.github.com/hydra-genetics/cnv_sv/commit/d5fdce851f04adb90f60f8588b387e99be8024bb))
* update docs/softwares.md ([688c081](https://www.github.com/hydra-genetics/cnv_sv/commit/688c08145e10f4cbed7c608cf346dc52353ec07b))
* update docs/softwares.md ([d76e841](https://www.github.com/hydra-genetics/cnv_sv/commit/d76e8411c74cec9f73e350bd48d9eb3254050720))
* update docs/softwares.md ([4268631](https://www.github.com/hydra-genetics/cnv_sv/commit/4268631f82c8fef281c584ccc25a9949f22f4ade))
* update the module output files ([02aca91](https://www.github.com/hydra-genetics/cnv_sv/commit/02aca91417bfabed6eb6ee1dc8583c609b9cc02b))

## [0.9.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.8.0...v0.9.0) (2025-05-14)


### Features

* added scanITD ([6b7e9bf](https://www.github.com/hydra-genetics/cnv_sv/commit/6b7e9bf5a6ecf48c518f4542f9a6bfb2db383661))


### Bug Fixes

* optional bed file ([7e2d882](https://www.github.com/hydra-genetics/cnv_sv/commit/7e2d8821e7206859d58532f0889455f2abf01ef3))
* **scanITD:** correct rule params ([01cb157](https://www.github.com/hydra-genetics/cnv_sv/commit/01cb157b762dd03023eea58aefcadd74a2fde9a6))
* **scanITD:** correct shell command ([3da1a76](https://www.github.com/hydra-genetics/cnv_sv/commit/3da1a76bf95766149f9d750f9472d8ec3ed1b185))
* update to bugfixed docker ([20477be](https://www.github.com/hydra-genetics/cnv_sv/commit/20477be5d70c92803b258398089e478015289f6b))


### Documentation

* added docs and tests ([0fc6464](https://www.github.com/hydra-genetics/cnv_sv/commit/0fc64646e1049d93a85021ce1a890220d37582f9))

## [0.8.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.7.1...v0.8.0) (2025-05-05)


### Features

* add entries for severus_t_only and severus_tn ([89ec119](https://www.github.com/hydra-genetics/cnv_sv/commit/89ec119182c830a3280a5db585abd80468ebc883))
* add function to compile input BAM paths ([49c2f1c](https://www.github.com/hydra-genetics/cnv_sv/commit/49c2f1c15452863db3fa9990a1cff8f7d98491ca))
* add paraphase ([2d533ab](https://www.github.com/hydra-genetics/cnv_sv/commit/2d533ab7eda83cc077803e9751f33ccae60ac05c))
* add pbsv.smk ([3854216](https://www.github.com/hydra-genetics/cnv_sv/commit/385421622d3a7a7590b9b33fbc80515a55f3e322))
* add rule for HiFiCNV ([930b800](https://www.github.com/hydra-genetics/cnv_sv/commit/930b8008ea60cc36a3d9fe45d34018909aa08442))
* add sawfish ([e7408b8](https://www.github.com/hydra-genetics/cnv_sv/commit/e7408b8d30c29613e389123c5fe4c87dbfa9637a))
* add severus related entries to config.yaml ([cd64c98](https://www.github.com/hydra-genetics/cnv_sv/commit/cd64c9802e938117f41e2233a31a363164611f31))
* add severus.smk ([53c64dc](https://www.github.com/hydra-genetics/cnv_sv/commit/53c64dc0924901482bd85cb1f6f3b39dbc6b9137))
* add severus.smk to the list of included rules ([39a0dd0](https://www.github.com/hydra-genetics/cnv_sv/commit/39a0dd0208f1f9bf33b3b99a899880e9a091e869))
* include hificnv.smk and pbsv.smk ([76ed3fc](https://www.github.com/hydra-genetics/cnv_sv/commit/76ed3fc59f46b5cc661b97c06f1e67300bcd2203))
* new name and functionality of the functon for creating input bam files ([2b21dec](https://www.github.com/hydra-genetics/cnv_sv/commit/2b21dece2f6c8c38d3d336258c4a70793327e243))
* update severus version to 1.5 ([e227d4f](https://www.github.com/hydra-genetics/cnv_sv/commit/e227d4f37155ecf6f7e3fb09da472b3b92f16659))
* use get_input_bam() in rule cnvkit_batch ([db571ba](https://www.github.com/hydra-genetics/cnv_sv/commit/db571ba3ed66f13b057e0abf575dd2afd89e8101))


### Bug Fixes

* add 'exclude' entry and 'exclude' file for hificnv ([f9c9bf5](https://www.github.com/hydra-genetics/cnv_sv/commit/f9c9bf50cba19fc70a08bf4ae5c39e4a7d452a45))
* add {type} to input, output, log and benchmark files in pbsv_call ([6c5c78c](https://www.github.com/hydra-genetics/cnv_sv/commit/6c5c78cf1b7417988114d65e8b24db7e9dad6944))
* add {type} to ouput, log and benchmark files in pbsv_discover ([f4ca49c](https://www.github.com/hydra-genetics/cnv_sv/commit/f4ca49c69802b0eb7656e1ffa90eab66eba417de))
* add {type} to output, log and benchmark files in rule hificnv ([26b5c38](https://www.github.com/hydra-genetics/cnv_sv/commit/26b5c38a4aac977af85a3b335389fabb95d75abf))
* add {unit_type} to output files of pbsv_call, pbsv_discover and hificnv in compile_output_list() ([71b2572](https://www.github.com/hydra-genetics/cnv_sv/commit/71b2572f902ae8db0a54cb87b5eb1b606dcb825d))
* add all output files description to severus_* ([8705dce](https://www.github.com/hydra-genetics/cnv_sv/commit/8705dce34522631f47433085744f25be922e8428))
* add all required sub-entries to severus entries ([ee36340](https://www.github.com/hydra-genetics/cnv_sv/commit/ee36340ef96d5fcce60590c001824e86b17d82d8))
* add container for hificnv to .tests/config_pacbio.yaml ([c1fa21a](https://www.github.com/hydra-genetics/cnv_sv/commit/c1fa21a1b4633dc4970921fb474954a649d7ec47))
* add docker image path, reference and trf files to pbsv rules ([4842a93](https://www.github.com/hydra-genetics/cnv_sv/commit/4842a93e128db7f723ebea5a0138eceea81a198f))
* add entries for hificnv and pbsv ([624d28b](https://www.github.com/hydra-genetics/cnv_sv/commit/624d28b16292af848c290d85364c8d4b54e649d2))
* add output files from hificnv and pbsv to compile_output_list() ([d7241dd](https://www.github.com/hydra-genetics/cnv_sv/commit/d7241dd6c3a2199b07d561504e97d4f2f89bba87))
* add renaming of the output files after HiFiCNV run ([762fbd7](https://www.github.com/hydra-genetics/cnv_sv/commit/762fbd79d0f9ecc74377c585cac505e657d97d62))
* add required quotation marks in awk command ([ccdc1fd](https://www.github.com/hydra-genetics/cnv_sv/commit/ccdc1fd373c3bb5ac3b709b0e73582d61f3f366d))
* add temp() directive to output ([0c4add0](https://www.github.com/hydra-genetics/cnv_sv/commit/0c4add03f30a1322a1efdc41a29e3dc585aa9774))
* change input BAM files in pbsv_discover from hard-coded to function ([cab49eb](https://www.github.com/hydra-genetics/cnv_sv/commit/cab49ebb15aecac2e6b53d62b4af7b333b54f00c))
* change input bam path; mv ref to params ([6075a12](https://www.github.com/hydra-genetics/cnv_sv/commit/6075a12fdfd8d355c1ee391b8168616ec9d9fe7f))
* change input from  path to function ([8b7c139](https://www.github.com/hydra-genetics/cnv_sv/commit/8b7c139bff70774c2953654ba559c72d7a16bb6f))
* change output/log to vcf.gz ([7d409db](https://www.github.com/hydra-genetics/cnv_sv/commit/7d409dbf66e7e74cf6204e8e3a14afa9623dd47a))
* correct alignment_path and index_path in get_input_bam() ([d232314](https://www.github.com/hydra-genetics/cnv_sv/commit/d2323143815b5bc19854817f8bbd4fdb4a539b38))
* correct config entry name for pon ([c867311](https://www.github.com/hydra-genetics/cnv_sv/commit/c867311cd4e644b78ea21c10386bb5010ea3f25b))
* correct output subdir name for severus_tn ([74c0648](https://www.github.com/hydra-genetics/cnv_sv/commit/74c064801b0ab0efcb7c1a7b10e476852a36aeb0))
* Delete workflow/rules/juli.smk ([d72fc1e](https://www.github.com/hydra-genetics/cnv_sv/commit/d72fc1e5dd68333144439d86cb24fe8c4abc1abe))
* limit number of threads used by Jumble to that specified in config ([b27fb17](https://www.github.com/hydra-genetics/cnv_sv/commit/b27fb174baa146d5492e787232f6cb4a3805477f))
* make two entries for pbsv (call & discover) ([c8a9272](https://www.github.com/hydra-genetics/cnv_sv/commit/c8a9272e129f1bd49c9b355ada28353dc075995f))
* make two entries for pbsv (call & discover) ([55b8793](https://www.github.com/hydra-genetics/cnv_sv/commit/55b8793b11a3b0688a42beea99474c458251bdfa))
* **paraphase:** add missing extra to shell command ([1250f24](https://www.github.com/hydra-genetics/cnv_sv/commit/1250f243de2623e1c235a26258ac32febceafa45))
* remove print() from get_input_bam() ([80800a4](https://www.github.com/hydra-genetics/cnv_sv/commit/80800a4e1724d0920f6a6ea7a1487f2ec3fa57f7))
* **sawfish:** update inputs to function ([671b316](https://www.github.com/hydra-genetics/cnv_sv/commit/671b316e8f2c7410c08447f0ae959e3862bbd9c0))
* separate output files not only dir ([f977330](https://www.github.com/hydra-genetics/cnv_sv/commit/f97733055ef615689608e95680864fa528b69058))
* **sniffles2:** add sample-id param to add sample name in VCF ([7f236b9](https://www.github.com/hydra-genetics/cnv_sv/commit/7f236b990664d66de412663bc4cf4913836187bf))
* specify all output files and fix the 'ouput file name as prefix' problem ([8d9490e](https://www.github.com/hydra-genetics/cnv_sv/commit/8d9490ef4c232a1c0dcf8124b22f0af79a612ae0))
* **tiddit:** add threads and access extras ([330c506](https://www.github.com/hydra-genetics/cnv_sv/commit/330c5064b76edbb8205b6d630f311a3dae3ff4ca))
* use correct name of the input bam function ([0ff7ccd](https://www.github.com/hydra-genetics/cnv_sv/commit/0ff7ccdf1c4eb70f3990fb25bf960c848e2c7de6))
* use get_input_bam() ([5aeb9af](https://www.github.com/hydra-genetics/cnv_sv/commit/5aeb9af6deb282743a7191430f78bd8b63ccf85a))


### Documentation

* update schemas with sawfish info ([e8f0163](https://www.github.com/hydra-genetics/cnv_sv/commit/e8f0163b632499f4f81908dc5d76d7d24c48e763))

### [0.7.1](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.7.0...v0.7.1) (2024-10-16)


### Bug Fixes

* add missing parenthesis ([9328c59](https://www.github.com/hydra-genetics/cnv_sv/commit/9328c593d2bc0b3b525ef1fd8214b90e9ed2ee4a))
* jumble results removed prematurely ([5ae3b74](https://www.github.com/hydra-genetics/cnv_sv/commit/5ae3b74b4da103c03332c77f2caa19aed5e2bee1))
* **tiddit:** add bai file as an input ([ace4bcd](https://www.github.com/hydra-genetics/cnv_sv/commit/ace4bcd17b2f7c16f58b9eac6b223308067a9d45))


### Documentation

* update documentation ([087f67e](https://www.github.com/hydra-genetics/cnv_sv/commit/087f67e52ce4cda10861f58f4c42609b354dad95))
* update mkdocs ([e28d589](https://www.github.com/hydra-genetics/cnv_sv/commit/e28d5896c9dbd3727bcff704bf6e7967b2b43bf0))

## [0.7.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.6.0...v0.7.0) (2024-10-02)


### Features

* add possibility to change caller annotation ([267285f](https://www.github.com/hydra-genetics/cnv_sv/commit/267285f782855a7903c4165ef829c515d2082aa2))
* add possibility to change caller name ([0b5f000](https://www.github.com/hydra-genetics/cnv_sv/commit/0b5f000f761b934d4ffc747b432f00efdffdebec))
* make caller name configurable ([461aee1](https://www.github.com/hydra-genetics/cnv_sv/commit/461aee1d4d14bc2cf337c8c634eebe2148827325))


### Bug Fixes

* make name configurable ([69ea10a](https://www.github.com/hydra-genetics/cnv_sv/commit/69ea10ac806fd0fec770a72e3974d7e1b6a18a14))
* put param in correct rule ([b99712e](https://www.github.com/hydra-genetics/cnv_sv/commit/b99712e44513dc7840a1da9fac5b0244ecb60b15))
* set sample_name correctly in header ([ef428c6](https://www.github.com/hydra-genetics/cnv_sv/commit/ef428c6d69394fd8495eaf8ebfe5db730932f475))

## [0.6.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.5.0...v0.6.0) (2024-09-12)


### Features

* add chr to chromosome if missing ([0197bb2](https://www.github.com/hydra-genetics/cnv_sv/commit/0197bb2c48b9d45d72c9d562a18ac02169511508))
* add sniffles2 single sample calling ([0ce9b8a](https://www.github.com/hydra-genetics/cnv_sv/commit/0ce9b8aa3821e4c98d78d8f89ca8d423164e1765))
* add tabix and bgzip for vcf output ([f29e4a2](https://www.github.com/hydra-genetics/cnv_sv/commit/f29e4a29c469b931137f64a9698ccdb3c786c744))
* add trgt ([145641d](https://www.github.com/hydra-genetics/cnv_sv/commit/145641daefdd6e47d9dc07b229c964582f8874d6))
* added cnvkit call of jumble output ([605ce90](https://www.github.com/hydra-genetics/cnv_sv/commit/605ce9040f4137aead0a1f392d34cc0d931ee5fe))
* added jumble_vcf rule for segment exports ([0aa2632](https://www.github.com/hydra-genetics/cnv_sv/commit/0aa263211e11c6f16c56b1adb47da832cb87b905))
* added the CNV caller Jumble ([4d24aa8](https://www.github.com/hydra-genetics/cnv_sv/commit/4d24aa883305ac8783af0df91f9b9e02e2071545))


### Bug Fixes

* add default values for pacbio_alignment and ont_alignment in config.schema.yaml ([ddedade](https://www.github.com/hydra-genetics/cnv_sv/commit/ddedade1728750d9cdd64f18a4d53532a05457ac))
* change to compatible wrapper ([4759cff](https://www.github.com/hydra-genetics/cnv_sv/commit/4759cffbb9cbd2ac24efb85296cb9a854fcc19ef))
* code improvement ([911422a](https://www.github.com/hydra-genetics/cnv_sv/commit/911422a2cd2428c574b7f398c8129f0a340577a6))
* **common:** set type to str for pandas dataframe ([4c0e0f2](https://www.github.com/hydra-genetics/cnv_sv/commit/4c0e0f2252da1ab187ff94c132523fe1958c959c))
* correct input parameter name ([1722db4](https://www.github.com/hydra-genetics/cnv_sv/commit/1722db4fea07392d23eb31c13bf6bd4586fb5ae6))
* corrected output file names ([238fb36](https://www.github.com/hydra-genetics/cnv_sv/commit/238fb36ac6e9ffe4bb6318ef244e564e5f7719ce))
* jumble_vcf should have same parameters as cnvkit_vcf ([1748415](https://www.github.com/hydra-genetics/cnv_sv/commit/1748415286965c99ea7864c5fc11694f3f8f6568))
* make sure sample column is a string ([2d3aa6b](https://www.github.com/hydra-genetics/cnv_sv/commit/2d3aa6bae5addca7219e9aa833778110efd1173f))
* missing parenthesis ([be912f8](https://www.github.com/hydra-genetics/cnv_sv/commit/be912f822c7f7284676c7b083dc51214da5bb9d8))
* update common container ([409dbda](https://www.github.com/hydra-genetics/cnv_sv/commit/409dbda6607a2d6cbba60ea63632bebb9893e333))
* Update config_pacbio.yaml ([e7f6117](https://www.github.com/hydra-genetics/cnv_sv/commit/e7f6117eac6d519f1238395b0bb79864b6d9daa9))
* Update config.yaml ([c2a82f8](https://www.github.com/hydra-genetics/cnv_sv/commit/c2a82f81e283dabec54a7e8494ba1adf41a2b574))
* Update config.yaml ([ea9a0f7](https://www.github.com/hydra-genetics/cnv_sv/commit/ea9a0f7d916952197381fcb61c98063a2f71835b))
* Update config.yaml ([794e3d6](https://www.github.com/hydra-genetics/cnv_sv/commit/794e3d6fedc8d26674128bb7f06db814e7900d44))
* wrong config parameter name ([3ad5336](https://www.github.com/hydra-genetics/cnv_sv/commit/3ad533667c4ea0577793c7e55cdb3cf3474cb23e))


### Documentation

* documentation fix ([958550b](https://www.github.com/hydra-genetics/cnv_sv/commit/958550be4c0b91ac499f73310ce4f4078e76b6b0))
* update rules ([e1166af](https://www.github.com/hydra-genetics/cnv_sv/commit/e1166aff1bee8ace4b6f4a8583cfd9d7a4cd7218))

## [0.5.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.4.1...v0.5.0) (2024-03-28)


### Features

* add names for the vcfs using params and add bnd_distance param ([e32ca47](https://www.github.com/hydra-genetics/cnv_sv/commit/e32ca47bfbac69ef37bf670c9a014937c0cfd037))
* add priority as a param ([db6e1ba](https://www.github.com/hydra-genetics/cnv_sv/commit/db6e1bac37321f23a7d2763b94493ebc0f4dd049))
* make the tumor cell content optional ([#160](https://www.github.com/hydra-genetics/cnv_sv/issues/160)) ([1c960c6](https://www.github.com/hydra-genetics/cnv_sv/commit/1c960c638ffdefea8b23296f8a49fbe65b82b557))


### Bug Fixes

* add missing space in command ([21e2767](https://www.github.com/hydra-genetics/cnv_sv/commit/21e2767e8741a6d06cb2b4fedc5928d21025f49c))
* add quite to picard for stdout ([8f09906](https://www.github.com/hydra-genetics/cnv_sv/commit/8f09906109194e6444e84c7b3fc6a41ecfebab73))
* added RTD for all reviewer rules and fix rule name ([3a81357](https://www.github.com/hydra-genetics/cnv_sv/commit/3a81357f9ce3d063d9a77403e06d1ea6a557b26a))
* allow missing TC when generating GATK VCF ([#162](https://www.github.com/hydra-genetics/cnv_sv/issues/162)) ([aa6fece](https://www.github.com/hydra-genetics/cnv_sv/commit/aa6fece935ea8587910b273da7cbf6b8674bfb68))
* bump snakemake version to support apptainer ([#151](https://www.github.com/hydra-genetics/cnv_sv/issues/151)) ([eb4b05d](https://www.github.com/hydra-genetics/cnv_sv/commit/eb4b05df754c2b6be16d5feed3244d00d9f4b5f1))
* **cnvpytor:** add bam index as input ([60c840a](https://www.github.com/hydra-genetics/cnv_sv/commit/60c840ad7c71dbfd66c0925becc210ae31dfa41a))
* handle that cns files can have different number of columns ([57da77b](https://www.github.com/hydra-genetics/cnv_sv/commit/57da77b993ea25be88116e75712b09411e45677e))
* incorrect type for copy number thresholds ([#164](https://www.github.com/hydra-genetics/cnv_sv/issues/164)) ([b549266](https://www.github.com/hydra-genetics/cnv_sv/commit/b549266994149bf7b4bdcc385456229b4f0eee77))
* pin pulp version to <2.8.0 ([1c960c6](https://www.github.com/hydra-genetics/cnv_sv/commit/1c960c638ffdefea8b23296f8a49fbe65b82b557))
* rm extra } in benchmark name ([e2f34ce](https://www.github.com/hydra-genetics/cnv_sv/commit/e2f34ce8d24515d7d1f980994f3d72e12878d4da))
* specify model when calling and filtering in cnvpytor, and threads in readdpth rule ([4347806](https://www.github.com/hydra-genetics/cnv_sv/commit/434780660b693fc7b9dbdf709cc4700be33de63c))
* **svdb:** use the overlap and extra params ([fe656a9](https://www.github.com/hydra-genetics/cnv_sv/commit/fe656a91dd0a74c566ef0b8bf3512609a071ab94))
* update container for rule pindel_update_vcf ([#161](https://www.github.com/hydra-genetics/cnv_sv/issues/161)) ([512a924](https://www.github.com/hydra-genetics/cnv_sv/commit/512a924beba54282757745c8c66475e78f52cca3))


### Performance Improvements

* add bcftools reheader to pindel vcf for correct samplename ([fcadaba](https://www.github.com/hydra-genetics/cnv_sv/commit/fcadaba56976bad8a6f2975558c5368c25016fb8))


### Documentation

* added RTD files and the rules for automap and cnvkit ([ce2f2ca](https://www.github.com/hydra-genetics/cnv_sv/commit/ce2f2ca0b4cdd6876968f638d3d30226f62bfe19))
* added RTD for all manta rules ([2571744](https://www.github.com/hydra-genetics/cnv_sv/commit/25717444f3a6a48825e95d2c4e756ec2af62c244))
* added RTD for all pindel rules ([b0ca4da](https://www.github.com/hydra-genetics/cnv_sv/commit/b0ca4dadbbada7bbad66c1274c125f58ef69447d))
* added RTD for all purecn rules ([7211b29](https://www.github.com/hydra-genetics/cnv_sv/commit/7211b298c138c8aaa693bca072b78885379d8532))
* added RTD for all smncaller rules ([ed55b6b](https://www.github.com/hydra-genetics/cnv_sv/commit/ed55b6b121c08fb4b4935fe231f493e563f41104))
* added RTD for all svdb rules ([1cdd6ba](https://www.github.com/hydra-genetics/cnv_sv/commit/1cdd6ba37ca9cfd3c21458188d3d04b706334112))
* added RTD for cnvpytor and exomedepth ([8fd490a](https://www.github.com/hydra-genetics/cnv_sv/commit/8fd490a8b73efc412797859678569ecba07ac783))
* added RTD for expansionhunter ([4f290ca](https://www.github.com/hydra-genetics/cnv_sv/commit/4f290ca1e0b69dc7540411d877d380644da69b7e))
* added RTD for gatk cnv rules ([a0a575e](https://www.github.com/hydra-genetics/cnv_sv/commit/a0a575ef4a37ac0014b67da8690f8241eeef68f9))
* added RTD for tiddit rules ([fea5964](https://www.github.com/hydra-genetics/cnv_sv/commit/fea5964ec64f5e184c3df75a43241485fbdfa8b0))
* added RTD for upd rules ([2eccfe5](https://www.github.com/hydra-genetics/cnv_sv/commit/2eccfe5bc01c485afc3996cd85cf0c3700a06175))
* dag graph ([7a654c2](https://www.github.com/hydra-genetics/cnv_sv/commit/7a654c2ddedecc9a47321f8120b955162fe66674))
* organise rules based on software ([36568d9](https://www.github.com/hydra-genetics/cnv_sv/commit/36568d9cc61789b3ee5e07af1d695c01d783e496))
* Update docs/softwares.md ([d3e3f67](https://www.github.com/hydra-genetics/cnv_sv/commit/d3e3f67f1e91e15223cc655133c944128c9e5a54))
* Update docs/softwares.md ([188ef59](https://www.github.com/hydra-genetics/cnv_sv/commit/188ef594302900cc0b37580e616775412fe79764))
* Update docs/softwares.md ([90ad983](https://www.github.com/hydra-genetics/cnv_sv/commit/90ad983f53a88117f0532bb8042c74c279d1d7d5))
* Update docs/softwares.md ([dd02db5](https://www.github.com/hydra-genetics/cnv_sv/commit/dd02db57a7d41c253c1e1d5e79b521d43a6f003d))
* Update docs/softwares.md ([6f16c13](https://www.github.com/hydra-genetics/cnv_sv/commit/6f16c13479c785cfa52f7c318a9c85aa6a2f5727))
* Update docs/softwares.md ([25b3195](https://www.github.com/hydra-genetics/cnv_sv/commit/25b31959177b1cd798e4896e47c5bab36c94bad3))
* Update docs/softwares.md ([c02990c](https://www.github.com/hydra-genetics/cnv_sv/commit/c02990c978a76a7d26ed8f1ce4f962e42bc566e3))
* update to new rule plugin ([aedbdfd](https://www.github.com/hydra-genetics/cnv_sv/commit/aedbdfd8d9f3084b4e49a0d7e9bb709da707fb10))
* Update workflow/schemas/resources.schema.yaml ([6136f71](https://www.github.com/hydra-genetics/cnv_sv/commit/6136f71519cb149297621d16cb16c31bf8017838))

### [0.4.1](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.4.0...v0.4.1) (2023-05-04)


### Bug Fixes

* change to config.get in exomedepth_call ([34f5753](https://www.github.com/hydra-genetics/cnv_sv/commit/34f5753530450e6ae8ba51bd19c5bb0946a45d14))

## [0.4.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.3.1...v0.4.0) (2023-05-02)


### Features

* add upd ([cc5a660](https://www.github.com/hydra-genetics/cnv_sv/commit/cc5a660a61d0f6d9b4a8a4bc6559cadfcf26f741))
* **exomedepth:** add checkpoint rule to allow reference file as input ([a75d1e5](https://www.github.com/hydra-genetics/cnv_sv/commit/a75d1e5939e994ba44e665d8f6ae537228e71b5c))
* **exomedepth:** Add option to run exomedepth with hg38 ([179bab3](https://www.github.com/hydra-genetics/cnv_sv/commit/179bab3da18e62a927b6b93d2daef7727ebdd0b6))
* remove conda support and testing ([0a3610d](https://www.github.com/hydra-genetics/cnv_sv/commit/0a3610d881bd6c5a95dc13754140220b6adf3c13))


### Bug Fixes

* automap input and integration test ([20ec062](https://www.github.com/hydra-genetics/cnv_sv/commit/20ec062be15dce5a3766f50254d3fdd73b86361e))
* **exomedepth:** chosse female ref when the peddy sex is NA ([4357702](https://www.github.com/hydra-genetics/cnv_sv/commit/4357702e36f4b472b57bf6643ea19f26bf6356b0))
* get output ([fb11ed2](https://www.github.com/hydra-genetics/cnv_sv/commit/fb11ed2415ae8e281adb6000b531eb7974ed1bea))
* get samples filename from config, close [#135](https://www.github.com/hydra-genetics/cnv_sv/issues/135) ([46b53c6](https://www.github.com/hydra-genetics/cnv_sv/commit/46b53c6f51a040ff12e2a3f6c74ff55a1e5379ff))
* linting ([e28d277](https://www.github.com/hydra-genetics/cnv_sv/commit/e28d27753c0373190050474ecee1b5637151f596))
* out_dir ([55a13f0](https://www.github.com/hydra-genetics/cnv_sv/commit/55a13f040a0ccf8a8ab24e02e3af3efa9e08bb82))
* output ([8b8af48](https://www.github.com/hydra-genetics/cnv_sv/commit/8b8af48e9e0b47224fb028eeff34850431ad7df9))
* tests for sutomap ([33a384e](https://www.github.com/hydra-genetics/cnv_sv/commit/33a384e5634d93a67cdb9b910ecfb632d40a38c1))
* tsv to pdf ([3aa56c7](https://www.github.com/hydra-genetics/cnv_sv/commit/3aa56c72aa6d2d626b08de057310cc264fec7235))
* typo ([d5cb0ca](https://www.github.com/hydra-genetics/cnv_sv/commit/d5cb0ca43f30e8b923a21dadc96c04b7c44b09f3))


### Documentation

* update compatibility ([0e1e381](https://www.github.com/hydra-genetics/cnv_sv/commit/0e1e381ac649bb48a0573f932053d48153e86beb))
* update compatibility ([fe31a4b](https://www.github.com/hydra-genetics/cnv_sv/commit/fe31a4b4b39cb5413c8875947d22156a416ea1e2))
* update disclaimer ([1d98fea](https://www.github.com/hydra-genetics/cnv_sv/commit/1d98feab3ab2f1c6be40909119c00baee9045a48))
* update disclaimer ([cbb7592](https://www.github.com/hydra-genetics/cnv_sv/commit/cbb7592ecf98ffb4245e7673bcbefde2a4759a82))
* update message format ([76eaa0b](https://www.github.com/hydra-genetics/cnv_sv/commit/76eaa0bca03ee6942ee83b63dd68f5bf256f364d))
* update message format ([4e4f73a](https://www.github.com/hydra-genetics/cnv_sv/commit/4e4f73ad6f3c3e2ba5cd0896475ee63e8ebe5533))
* updated output files ([d56c8d6](https://www.github.com/hydra-genetics/cnv_sv/commit/d56c8d6282382719dc2ce0c5751a2097f616f8b6))

### [0.3.1](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.3.0...v0.3.1) (2023-01-31)


### Bug Fixes

* **purecn_copy_output:** add missing benchmark and params ([c10d758](https://www.github.com/hydra-genetics/cnv_sv/commit/c10d7587ac96403dbfbee6df4cd3436a08377df3))


### Documentation

* update compatibility file ([b9fd87d](https://www.github.com/hydra-genetics/cnv_sv/commit/b9fd87d767499226c34dedc11aebf73fa95d2152))

## [0.3.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.2.0...v0.3.0) (2023-01-27)


### Features

* added default in function ([bf4f3f2](https://www.github.com/hydra-genetics/cnv_sv/commit/bf4f3f2399494742c9ed38ab24cd2fa619dde4d5))
* added function that extracts purity for purecn ([9cf298e](https://www.github.com/hydra-genetics/cnv_sv/commit/9cf298e8eda05dc02280a84aa440dafc5f59082f))
* added pathology option ([4b3adc5](https://www.github.com/hydra-genetics/cnv_sv/commit/4b3adc52718494c0673b1a8bf7252e9bcedc0295))
* added rule purecn_purity_file ([2bed973](https://www.github.com/hydra-genetics/cnv_sv/commit/2bed973dd8e9e2420d630978fe24571fec36a388))
* added tags for purity method ([b9f3d3d](https://www.github.com/hydra-genetics/cnv_sv/commit/b9f3d3d6b3a8ffbd39f32071d4eed12907703d6d))
* svdb is now configurable to tc method ([c3eb2e1](https://www.github.com/hydra-genetics/cnv_sv/commit/c3eb2e17e04eed00e8c0ea8314e59035560c4e4f))


### Bug Fixes

* added backslash in cmd ([994fa78](https://www.github.com/hydra-genetics/cnv_sv/commit/994fa783459192db0d6b039023a16a85907a4a13))
* allow purecn to run on all sample types ([6947acf](https://www.github.com/hydra-genetics/cnv_sv/commit/6947acf0a3e91aa3de0080cb9a907fe8af10b880))
* bugfx ([ff30237](https://www.github.com/hydra-genetics/cnv_sv/commit/ff30237dec88ebf163dffea7f8d55e930968c6e4))
* corrected get_tc function output ([3c89633](https://www.github.com/hydra-genetics/cnv_sv/commit/3c89633a26be3a315774b112640d467d5bb15705))
* input file now with germline annotation ([bb7e5dd](https://www.github.com/hydra-genetics/cnv_sv/commit/bb7e5ddb7c1bbb0305d4df546cd8d43002a8bc21))
* output files ([97fb714](https://www.github.com/hydra-genetics/cnv_sv/commit/97fb714803894e91600badb34b9a89b4483992b0))
* rename caller from gatk_cnv to gatk ([6c0d00f](https://www.github.com/hydra-genetics/cnv_sv/commit/6c0d00f149419c883935b80faa30d17ca483ecce))
* rule name and output folder ([003f415](https://www.github.com/hydra-genetics/cnv_sv/commit/003f4159f3a3a681ef393a94cde2ec1b5012784c))
* tc files only in params and does not trigger purecn ([28bdbe4](https://www.github.com/hydra-genetics/cnv_sv/commit/28bdbe4e99c6569b20cd33eaf71ad6d9be00205d))


### Documentation

* update CODEOWNERS ([c20ccbc](https://www.github.com/hydra-genetics/cnv_sv/commit/c20ccbcdee442834e139f5e510f92425b7cdad68))
* update compatibility ([6bb06a5](https://www.github.com/hydra-genetics/cnv_sv/commit/6bb06a5624d813f0a67b06bc4642a007c2bc4ffa))

## [0.2.0](https://www.github.com/hydra-genetics/cnv_sv/compare/v0.1.0...v0.2.0) (2022-11-22)


### Features

* Initial PureCN implementation ([827ecd5](https://www.github.com/hydra-genetics/cnv_sv/commit/827ecd59ce1f38d50e75eae16f58dcd5fdc3e906))
* run expansionhunter with sex from peddy ([96248bf](https://www.github.com/hydra-genetics/cnv_sv/commit/96248bfb2dd5f6625905d53de20c87ed4be9955d))
* update pindel call wrapper ([deca8ee](https://www.github.com/hydra-genetics/cnv_sv/commit/deca8eedbf4effddaafe51964d9760f9245e18e2))


### Bug Fixes

* Add cnvkit_seg output ([6b2bc02](https://www.github.com/hydra-genetics/cnv_sv/commit/6b2bc024ebd445f0a44ae10c9e09b8f0b5c4682c))
* Add dummy intervals file for dry-run ([c36efa8](https://www.github.com/hydra-genetics/cnv_sv/commit/c36efa8e753df8b3875ad766b69099d400e41ca8))
* Add mutect2 VCFs for PureCN ([58c4784](https://www.github.com/hydra-genetics/cnv_sv/commit/58c4784727c27977ffa6683a198197da35caa68a))
* Add required parameters for purecn ([a6771a6](https://www.github.com/hydra-genetics/cnv_sv/commit/a6771a635b97209d17356e59e65d7d9182c56c82))
* fix integration test problem ([139c948](https://www.github.com/hydra-genetics/cnv_sv/commit/139c948a7a3ab963abcd6b4d83f914b35ba22fc6))
* fun_segmentation is a required property for purecn ([2376c99](https://www.github.com/hydra-genetics/cnv_sv/commit/2376c995fccc36dd3bae807ad12c098d60cd0877))
* make run step of manta require bam as input to prevent bam files from being cleaned away ([5918c73](https://www.github.com/hydra-genetics/cnv_sv/commit/5918c739006e1b7e4b25b11428ffd5e90801bad9))
* Minor language tweaks ([aebacab](https://www.github.com/hydra-genetics/cnv_sv/commit/aebacabd1f309d4d39699c298144e6dc28f3f3bd))
* Re-order resources in rules ([1c889d4](https://www.github.com/hydra-genetics/cnv_sv/commit/1c889d4bdbd71df099d8af5485ddbff8e74406e9))
* Rename cnvkit_seg log and benchmark ([c9f2808](https://www.github.com/hydra-genetics/cnv_sv/commit/c9f2808ee711a00a95ad766481e334d6cb3b2e13))
* Resolve snakefmt warnings ([fa7651e](https://www.github.com/hydra-genetics/cnv_sv/commit/fa7651e8b4d37db4b9e4d38ffea0cc1c75094773))
* Resolve snakefmt warnings ([b738559](https://www.github.com/hydra-genetics/cnv_sv/commit/b738559ec662e794f834b731b581afc633f41ea3))
* Separate PureCN output rule ([ef80f60](https://www.github.com/hydra-genetics/cnv_sv/commit/ef80f603a6d5068ad3ebb5c9fde82df7edc38b0e))
* Simplify expected PureCN output ([22fb008](https://www.github.com/hydra-genetics/cnv_sv/commit/22fb0089163a29416462ec60c2e3d7dce9dcbcb8))
* **tiddit:** fix env and bam to work with tiddit ([6803a15](https://www.github.com/hydra-genetics/cnv_sv/commit/6803a150934b149925919ce1f6bfdc19f49c83d6))
* Typos in output definitions ([106edaa](https://www.github.com/hydra-genetics/cnv_sv/commit/106edaabd90bd589140743aed519c175f4eff89c))
* Use dummy intervals ([8bcdae0](https://www.github.com/hydra-genetics/cnv_sv/commit/8bcdae0300b7ff849b5bd5015939d2c2595bf2dc))
* Use hydragenetics container for cnvkit_seg ([65bd30b](https://www.github.com/hydra-genetics/cnv_sv/commit/65bd30b030d8774190e488fdaaad5b37ebf17985))
* Use same conda env for all PureCN rules ([a795870](https://www.github.com/hydra-genetics/cnv_sv/commit/a795870b2c0628e9baeb2e331d13540376e43f36))


### Documentation

* **README:** change header to regular text ([d0e74d4](https://www.github.com/hydra-genetics/cnv_sv/commit/d0e74d4fa7dad72523c41d3a8209d40aceeaf0bf))

## 0.1.0 (2022-10-14)


### Features

* Add 3 rules for pindel plus envs and script ([5c51855](https://www.github.com/hydra-genetics/cnv_sv/commit/5c5185593b7718fa49270f017bd26e33f63b1243))
* add conventional-prs workflow ([94cc587](https://www.github.com/hydra-genetics/cnv_sv/commit/94cc5875c86cd2cdce2dc43a2a2cf4792337811c))
* add exomedepth ([ebfe037](https://www.github.com/hydra-genetics/cnv_sv/commit/ebfe037a87564692ee89be8e7114ddc959a6a4d3))
* Add logging and make script testable ([167a148](https://www.github.com/hydra-genetics/cnv_sv/commit/167a148bcf91823bbf0cde5f720bfb0c54c8fb7c))
* add missing entry in schema ([61a7055](https://www.github.com/hydra-genetics/cnv_sv/commit/61a70557f910f60ef3b9d7a4802f194c33f62cbb))
* Add rule cnvkit_call_loh which adapt calls according to TC ([3837d51](https://www.github.com/hydra-genetics/cnv_sv/commit/3837d51dfd94e11f4617fb7c5e887ea595376998))
* Add rule to generate cnvkit diagram ([bef8828](https://www.github.com/hydra-genetics/cnv_sv/commit/bef8828521497758e3c9f685dda009c6ca10a721))
* added bedfile to pindel call ([36aca30](https://www.github.com/hydra-genetics/cnv_sv/commit/36aca300ac068610384c0534bfe618b64d66aa5b))
* Added environments ([c919fed](https://www.github.com/hydra-genetics/cnv_sv/commit/c919fedec6e2b1601c8740fdbae8b178907a57ac))
* Added manta rules ([e77071d](https://www.github.com/hydra-genetics/cnv_sv/commit/e77071d8a702e85af9fcc73025bdc7bae29c7ba0))
* Added rule cnvkit coverage ([7162e00](https://www.github.com/hydra-genetics/cnv_sv/commit/7162e0084526b6782f986087f7ca5d735857a062))
* Added rule cnvkit_call together with test files ([4d7c808](https://www.github.com/hydra-genetics/cnv_sv/commit/4d7c80838866c49cd9dd33ffe124cd09df57e334))
* Added rule cnvkit_create_access ([7d5a550](https://www.github.com/hydra-genetics/cnv_sv/commit/7d5a550606a2f8134bf347e7b7ebca31a35bd887))
* Added rule cnvkit_create_antitargets ([5fd5d12](https://www.github.com/hydra-genetics/cnv_sv/commit/5fd5d12ea96a4f220823b93e233fa0ed9904ca4a))
* Added rule cnvkit_create_targets ([cdb3f5e](https://www.github.com/hydra-genetics/cnv_sv/commit/cdb3f5e14cd0de28173e48a2ebff8ff2a0ef3d8b))
* Added rule for manta tumor only calling. Refactor of old rule ([530a295](https://www.github.com/hydra-genetics/cnv_sv/commit/530a29503eecabc46fc814c7eb031fc0db6bb467))
* Added rule GATK_cnv_collectAllelicCounts ([d4abf5e](https://www.github.com/hydra-genetics/cnv_sv/commit/d4abf5ed996794397ee24cd3abd0640549e54d3b))
* Added rule GATK_cnv_collectReadCounts ([66ab8f7](https://www.github.com/hydra-genetics/cnv_sv/commit/66ab8f7e1163723f0fee38bab1c7ccc8d1befda4))
* Added rule GATK_cnv_denoiseReadCounts ([99822b5](https://www.github.com/hydra-genetics/cnv_sv/commit/99822b57882a5b17bdadae43d31445f279598371))
* Added rule GATK_cnv_modelSegments ([8c88d06](https://www.github.com/hydra-genetics/cnv_sv/commit/8c88d06e8d210117f2e38e6f7759d94332a1a9fc))
* Added rule germline_vcf ([9bbabc4](https://www.github.com/hydra-genetics/cnv_sv/commit/9bbabc42345a7cf416c2fc292287c213987611b8))
* Changed cnvkit_vcf to own script ([db9a85d](https://www.github.com/hydra-genetics/cnv_sv/commit/db9a85d03ba086ac4650f163d935e9640369abb8))
* Generate output list programmatically ([ef95803](https://www.github.com/hydra-genetics/cnv_sv/commit/ef9580374e6b9b72549183b37f8d3742f8b68e94))
* germline vcf now performed by filter vcf using pipeline config and override in Snakemake ([0a427b6](https://www.github.com/hydra-genetics/cnv_sv/commit/0a427b605f4ce1eea3fa8d8b0f1ff4dfcf3b8b80))
* make config.yaml location more flexible ([f6f4416](https://www.github.com/hydra-genetics/cnv_sv/commit/f6f4416ec0594aaa92e681266fa33085a095127d))
* make configfile/confgilefiles argument mandatory ([02b4f25](https://www.github.com/hydra-genetics/cnv_sv/commit/02b4f25baa9dfd27cfcbbc19a6b6fe293c1962ad))
* make module compatible with latest version of snv_indels and alignment module. ([8ea6383](https://www.github.com/hydra-genetics/cnv_sv/commit/8ea638395a2fbed2d8814747b748ca4766bbf4e6))
* make output more configurable ([dbd9bbf](https://www.github.com/hydra-genetics/cnv_sv/commit/dbd9bbf3b32aaa2f22f265dc5df22fffee5e214a))
* Make unnecessary manta output vcfs temp ([46688c2](https://www.github.com/hydra-genetics/cnv_sv/commit/46688c2a213415323142d439a7908acc3fa67390))
* **manta:** added extra params to manta ([e99460f](https://www.github.com/hydra-genetics/cnv_sv/commit/e99460f224c86057478ffb1f8fcdda30b5577926))
* new rule cnvkit_vcf that exports segment files to vcf ([9769357](https://www.github.com/hydra-genetics/cnv_sv/commit/976935749710ddcadb0da0385f02351716f98eac))
* New rule svdb ([e7349b4](https://www.github.com/hydra-genetics/cnv_sv/commit/e7349b4417a833876576a32a1fc1030257e9487e))
* Prepare pindel config script for unit testing ([13b0264](https://www.github.com/hydra-genetics/cnv_sv/commit/13b02646689787c6c34151e3822a0ae8749cbc8f))
* setup integration test for exomedepth ([7a30fc0](https://www.github.com/hydra-genetics/cnv_sv/commit/7a30fc0697bc401aa5e63b5f12e733d01a835d77))
* update snakemake-version ([830124c](https://www.github.com/hydra-genetics/cnv_sv/commit/830124cb81318ec9eb631806b4e82b5e0ea5af4b))
* update svdb version ([d6745f6](https://www.github.com/hydra-genetics/cnv_sv/commit/d6745f65ae7b84182af5a5f5e4aeca12f4590392))
* Upgrade pindel to wrapper v85.0.1 ([9048434](https://www.github.com/hydra-genetics/cnv_sv/commit/9048434f3c08f9db153e26e67702d3c63bacb2d5))


### Bug Fixes

* add containers and change from bcftools view to zcat in rule germinline_vcf ([5fa941f](https://www.github.com/hydra-genetics/cnv_sv/commit/5fa941f7aa8f34492a712bece844e725b727d177))
* add contigs to pindel-vcf ([4a52b00](https://www.github.com/hydra-genetics/cnv_sv/commit/4a52b0038e18c2301e80106eec4017e5b4a6cf3f))
* add library to install for exomedepth ([63d19c2](https://www.github.com/hydra-genetics/cnv_sv/commit/63d19c221ff679f7dda47976dee62a965ca2fdec))
* Added resources, rule name fix, rm testfiles, rm wp1 in config ([f408a14](https://www.github.com/hydra-genetics/cnv_sv/commit/f408a14af5893ae2e034a29cfb3d7b75a62fa448))
* better vcf annotation for filtering ([8aa6da8](https://www.github.com/hydra-genetics/cnv_sv/commit/8aa6da8bed5b35d88c86cb45809c2ef92996d445))
* bugfix ([140842c](https://www.github.com/hydra-genetics/cnv_sv/commit/140842c587a64e7b25183c0b58ad67b903a850b0))
* change run column name to flowcell in units. ([a69af37](https://www.github.com/hydra-genetics/cnv_sv/commit/a69af37ad1759e8e155527655d056e399494e7b6))
* change TC to tumor_content ([8c11ee6](https://www.github.com/hydra-genetics/cnv_sv/commit/8c11ee695e48331bd9a829b4e6808129eea67457))
* change to point to master ([41e6d16](https://www.github.com/hydra-genetics/cnv_sv/commit/41e6d166dddcd5d0bccdbaabb7b42b81ca981079))
* change type from float to number ([1731573](https://www.github.com/hydra-genetics/cnv_sv/commit/1731573fb98df66c6eb497bbc5851a1b82632521))
* Changed name again ([ea530f1](https://www.github.com/hydra-genetics/cnv_sv/commit/ea530f1c59da79dc95b71721d907e5247ea33c23))
* **cnvkit:** change python version ([f315187](https://www.github.com/hydra-genetics/cnv_sv/commit/f31518771ec929884e222977b5d7e31a50e72374))
* **cnvkit:** change to correct conda env file ([be3d3ba](https://www.github.com/hydra-genetics/cnv_sv/commit/be3d3ba5d3a9bca28daa23031f93684723182608))
* Correct runWorkflow script path and add container to config ([057b151](https://www.github.com/hydra-genetics/cnv_sv/commit/057b1511bf66fe80acb0dc0ceba34233c208477c))
* don't install snakemake using mamba ([dc08143](https://www.github.com/hydra-genetics/cnv_sv/commit/dc08143c1c088d91d8b09a46c92029e846fd56fe))
* float in header ([4e30d39](https://www.github.com/hydra-genetics/cnv_sv/commit/4e30d39b5470061211cbc1577722a797381a0f03))
* float in header ([7c0c233](https://www.github.com/hydra-genetics/cnv_sv/commit/7c0c23333eed9a0f1d5779d150445c8add3c9d7d))
* handling of multiple databases ([978d218](https://www.github.com/hydra-genetics/cnv_sv/commit/978d2185aa2a0a0c9b282ee6099872577ceb9ba2))
* info ids to long ([541bd8d](https://www.github.com/hydra-genetics/cnv_sv/commit/541bd8d4e8df4566872e82ed7b475482900370d0))
* int to float in vcf headers ([331da1d](https://www.github.com/hydra-genetics/cnv_sv/commit/331da1d554ecddec9c5e5be2f328a097020a9308))
* lock singularity version to prevent bug with latest version ([bfe978e](https://www.github.com/hydra-genetics/cnv_sv/commit/bfe978ea5e29d00fd6835e23b090f5eb5020d877))
* reference under rule ([92fd035](https://www.github.com/hydra-genetics/cnv_sv/commit/92fd03538d20ae7043606fd3f8594f77f3ddcec4))
* remove unused dependency. ([81f7ed6](https://www.github.com/hydra-genetics/cnv_sv/commit/81f7ed66695ed3a472bad33c6b5646dee25eb27c))
* revert changes ([e74b497](https://www.github.com/hydra-genetics/cnv_sv/commit/e74b497594df70e11e881f95c169d958fc740432))
* rule name change ([f617b06](https://www.github.com/hydra-genetics/cnv_sv/commit/f617b0655fe332b359407f7c3ed357a68e222d5c))
* set correct env file ([c81d30d](https://www.github.com/hydra-genetics/cnv_sv/commit/c81d30de64d4cc291e32b270c3a8a763c10d7fa5))
* set strict mode for conda ([d070bda](https://www.github.com/hydra-genetics/cnv_sv/commit/d070bdab490f726339e5a03263e674dd67b58612))
* update r script to handle empty results ([6bd463b](https://www.github.com/hydra-genetics/cnv_sv/commit/6bd463b141a14ba16669f59b816ad49b7580824f))
* used corrected cn instead of cn ([dc59207](https://www.github.com/hydra-genetics/cnv_sv/commit/dc5920740d9cc54f12e73ca5390fb096f2d2f2a8))


### Documentation

* add and fix schema for picard_update_vcf_sequence_dictionary ([9bba1e9](https://www.github.com/hydra-genetics/cnv_sv/commit/9bba1e9cb72bf3b0c3e966006e61372b3bb6971a))
* add configfile option to example ([da7ff6b](https://www.github.com/hydra-genetics/cnv_sv/commit/da7ff6b239181bf9ac13cae94a68c7fe93045220))
* add description to all config- and resources-schemas ([d7103d4](https://www.github.com/hydra-genetics/cnv_sv/commit/d7103d40fbfa5c102dfcc7b27fcba50e4c29e2eb))
