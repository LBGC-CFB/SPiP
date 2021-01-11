# Changelog
All notable changes to this project will be documented in this file.

## Unrelease

## [1.1] - 2021-01-08
- Use RData files instead of samtools to get sequences
- Add two options to select gene or transcript to analyze

## [1.0] - 2020-05-15
- Change the size of the window to scn for cryptic/de novo creation
- Update getRefSeqDatabase.r script
- Update splicing probabilities
- Update output text of SPiP
- Change VCF example file
- Add policy to the SPiP script

https://github.com/raphaelleman/SPiP/commit/e3ce1878b4a8f1ea30e2ddab778fbe646f6f8cef<br>
https://github.com/raphaelleman/SPiP/commit/9e677cb4776f24f08738a87e251f908c9ca4f1e8<br>
https://github.com/raphaelleman/SPiP/commit/84b6510d10ac03f9a8bbbe8132766d1d2e2cc5de<br>
https://github.com/raphaelleman/SPiP/commit/b17569d98c7446ad24531e2d34eef3c3f59bc9e0<br>
https://github.com/raphaelleman/SPiP/commit/cf5b8b0b1a35e013d09b9075be80c1bee12c9b84<br>
https://github.com/raphaelleman/SPiP/commit/573a4e42353aad0d4107d2fd8f75eb4379c60b7e<br>
https://github.com/raphaelleman/SPiP/commit/e87fd77f2d4b3c2f3201b8f6393a4423c41aa8f5<br>

## [0.6.1] - 2020-04-24
- Check if samtools is executalbe
- Add function to read fa.fai file to get the chromosome nomenclature

https://github.com/raphaelleman/SPiP/commit/e1b9965742760db49cdf3680438b92927dc748ba<br>

## [0.6](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.6.r) - 2020-03-06
- Add exonic information
- Add function to get output in VCF format
- Change the parallel package to foreach and doParallel
- Remove function to get sequence by web API
- Add a MIT license
- Add ASCII logo

https://github.com/raphaelleman/SPiP/commit/879568f540711b4c6e7f5f015708be6b4d645ef8<br>
https://github.com/raphaelleman/SPiP/commit/dece37f570907412658cf929b4c0516ed9cac33f<br>
https://github.com/raphaelleman/SPiP/commit/e756ad3662c8f5bde624b25661e2ea5c15d02e81<br>
https://github.com/raphaelleman/SPiP/commit/c92b3664e368b3abe50bc55bd555d555277fc20f<br>
https://github.com/raphaelleman/SPiP/commit/51230f9a2945f9952eebff0c5d8c507338c86b05<br>
https://github.com/raphaelleman/SPiP/commit/03fb95b24249220dc58774d1633ab1b0c9517a68<br>
https://github.com/raphaelleman/SPiP/commit/9dd9c034ed8ed539aac6347c302721043d847cfb<br>
https://github.com/raphaelleman/SPiP/commit/293cdb290d77f4df0492bff987b0b84c03d73c73<br>
https://github.com/raphaelleman/SPiP/commit/901a4e65e506a84646b29fa82f929457c66d8de9<br>
https://github.com/raphaelleman/SPiP/commit/295d93c19e570d28bb87f81c6c2520374c9bd226<br>

## [0.5.1] - 2020-01-16
- Corrected error from the convertion transcriptomic coordinates to genomic coordinates

https://github.com/raphaelleman/SPiP/commit/405e0c83aee735a56b9ca0192f387eebbf1ada51<br>
https://github.com/raphaelleman/SPiP/commit/ae2d60d381bf4047a2c82c6e08579a712b017bb5<br>
https://github.com/raphaelleman/SPiP/commit/ebeb5bba11dc08e3bb0689af16bdcddd38809ee4<br>

## [0.5](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.5.r) - 2019-12-17
- Update the output of SPiP + add gene names
- Create a new logo
- Update ESR threshold to -1.10

https://github.com/raphaelleman/SPiP/commit/71a8c45f6479deae87856ed4af89f38625ec2ab2<br>
https://github.com/raphaelleman/SPiP/commit/548230ed44d2f293e6e40a5764dead2d87e82dd9<br>
https://github.com/raphaelleman/SPiP/commit/3947bcc3a0892fc0a0f2ece445cf702959c43d7b<br>
https://github.com/raphaelleman/SPiP/commit/f734cb05634f8b71a81897be403d91588676cd1e<br>

## [0.4.2] - 2019-09-27
- Update probabilites of splicing alteration

https://github.com/raphaelleman/SPiP/commit/031ee967eb3bb0bc375adbc60fa11f2f21a1e9ca<br>

## [0.4.1] - 2019-07-23
- Update VCF functions to take into account the meat-information of VCF files used as input

https://github.com/raphaelleman/SPiP/commit/688305a19251d8cb89687c6bf130f012831042d5<br>

## [0.4](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.4.r) - 2019-07-01
- Remove options to choice tools used by SPiP
- Update the output of SPiP

https://github.com/raphaelleman/SPiP/commit/b8b7a69066dcc5b5c3e1e13294c8a038200fb784<br>
https://github.com/raphaelleman/SPiP/commit/5a4fc118cb6a5d8cd1a601aae19490e9305b565c<br>
https://github.com/raphaelleman/SPiP/commit/7969f71def1f9d8683bfa16a072d0a435e175e3d<br>
https://github.com/raphaelleman/SPiP/commit/6853376703e38bbe891d55c48df7318d38623e4b<br>

## [0.3.1] - 2019-05-24
- Corrected convertion transcriptomic coordinates to genomic coordinates

https://github.com/raphaelleman/SPiP/commit/ab615730064116888311a1da16f0dd120440faad<br>

## [0.3](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.3.r) - 2019-05-16
- Add meta-header to the output
- Update probabilities of splicing
- Add the flag 'de novo' for splice site creation
- Add the PDF manual of SPiP
- Corrected getRefSeqDatabase.r to add gene names

https://github.com/raphaelleman/SPiP/commit/559ba236c773209019e5d8ce2eef0d18c3efdd79<br>
https://github.com/raphaelleman/SPiP/commit/702984e32be464fc58b51201c68d979c9f8b7927<br>
https://github.com/raphaelleman/SPiP/commit/e4d988fac1f939e69cb9218c72e048bc96c172b2<br>

## [0.2.4] - 2019-05-07
- Add getRefSeqDatabase.r srcipt permitting to reload and adapte RefSaq database for SPiP
- Add checking function to known if the DNA sequence permits to calculate scores
- Change the README file
- Add VCF file as an example from 1000 Genomes data

https://github.com/raphaelleman/SPiP/commit/f9bfb6513165a7fba8e36b47c38fc307ad9c7202<br>
https://github.com/raphaelleman/SPiP/commit/c913f0aa775d2070248c9bf5834f3d6ad266eca9<br>
https://github.com/raphaelleman/SPiP/commit/96a0ea5cf7d33161c8b4c7a1e731d4363298e485<br>
https://github.com/raphaelleman/SPiP/commit/b37529138f8d9dcef905ff741bb8905f81fec467<br>
https://github.com/raphaelleman/SPiP/commit/c76fc6d369f7670c137353c1a5e98b0b7930d79a<br>
https://github.com/raphaelleman/SPiP/commit/3276a621afd4ff06bd52909a7067aa9e1b5aae83<br>

## [0.2.3] - 2019-04-29
- Corrected VCF reading functions

https://github.com/raphaelleman/SPiP/commit/749d13742157b3f3c4c8abc43f0d201342a291d1<br>

## [0.2.2] - 2019-04-04
- Corrected parallel functions

https://github.com/raphaelleman/SPiP/commit/fabfd133fb70ba6513221b2752e5a6b49654d2f2<br>

## [0.2.1] - 2019-03-28
- Add parallel options

https://github.com/raphaelleman/SPiP/commit/bd19031d2ed2dab160177aea41bf6ea76f3c2176<br>

## [0.2](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.2.r) - 2019-03-01
- Git creation and push on https://github.com/raphaelleman/SPiP

https://github.com/raphaelleman/SPiP/commit/ad5bdcb94e5a798dc5f80ba9fcc46187e1418ab1<br>
https://github.com/raphaelleman/SPiP/commit/0c681edc71ae85eaacb0771c1bd14bd11753bdb6<br>
https://github.com/raphaelleman/SPiP/commit/98c7198959e66e47bb629040175842810144f015<br>
https://github.com/raphaelleman/SPiP/commit/7f9d821588ec1ab3403683f9d6b6fcbe599243b3<br>
https://github.com/raphaelleman/SPiP/commit/39e6aa5b9d9958adb52821c60494c93c533e3541<br>
https://github.com/raphaelleman/SPiP/commit/4e9d2aa09ff4500cd693b9e699432afdc1021e38<br>

## [0.1](https://github.com/raphaelleman/SPiP/blob/master/Releases/SPiPv0.1.r) - 2019-02-19
- First version online only for Windows OS
