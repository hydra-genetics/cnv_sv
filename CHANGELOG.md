# Changelog

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
