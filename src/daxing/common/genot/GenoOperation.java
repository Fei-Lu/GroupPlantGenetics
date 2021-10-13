package daxing.common.genot;

import pgl.infra.dna.snp.BiSNP;
import java.util.BitSet;

/**
 * Modify from GenotypeOperation, improved performance, easy for me to customize
 * Methods to merge, split, select genotype tables
 */
public class GenoOperation {

    /**
     * Merge the second genotype into the first genotype. Return a new one.
     * <p> Two genotypes should have the same number of taxa in the same order.
     * @param gt1
     * @param gt2
     * @return Return null if the merging is not successful
     */
    public static GenoGrid mergeGenotypesBySite(GenoGrid gt1, GenoGrid gt2) {
        if (gt1.getTaxaNumber() != gt2.getTaxaNumber()) return null;
        int snpCount = gt1.getSiteNumber()+gt2.getSiteNumber();
        BitSet[][] bArray = new BitSet[snpCount][3];
        int cnt = 0;
        for (int i = 0; i < gt1.getSiteNumber(); i++) {
            for (int j = 0; j < gt1.genoSite[0].length; j++) {
                bArray[cnt][j] = gt1.genoSite[i][j];
            }
            cnt++;
        }
        for (int i = 0; i < gt2.getSiteNumber(); i++) {
            for (int j = 0; j < gt2.genoSite[0].length; j++) {
                bArray[cnt][j] = gt2.genoSite[i][j];
            }
            cnt++;
        }
        BiSNP[] nsnps = new BiSNP[snpCount];
        cnt = 0;
        for (int i = 0; i < gt1.getSiteNumber(); i++) {
            nsnps[cnt] = gt1.snps[i];
            cnt++;
        }
        for (int i = 0; i < gt2.getSiteNumber(); i++) {
            nsnps[cnt] = gt2.snps[i];
            cnt++;
        }
        GenoGrid gt = new GenoGrid(bArray, GenoGrid.Direction.BySite, gt1.taxa, nsnps);
        return gt;
    }

    /**
     * Merge the second genotype into the first genotype. Return a new one.
     * <p> Two genotypes should have the same number of taxa in the same order.
     * @param gt1
     * @param gt2
     * @return Return null if the merging is not successful
     */
    public static GenoRows mergeGenotypesBySite(GenoRows gt1, GenoRows gt2) {
        if (gt1.getTaxaNumber() != gt2.getTaxaNumber()) return null;
        int snpCount = gt1.getSiteNumber()+gt2.getSiteNumber();
        SiteGeno[] geno = new SiteGeno[snpCount];
        int cnt = 0;
        for (int i = 0; i < gt1.getSiteNumber(); i++) {
            geno[cnt] = gt1.geno[i];
            cnt++;
        }
        for (int i = 0; i < gt2.getSiteNumber(); i++) {
            geno[cnt] = gt2.geno[i];
            cnt++;
        }
        GenoRows gt = new GenoRows(geno, gt1.taxa);
        return gt;
    }

    /**
     * Merge the second genotype into the first genotype. Return a new one.
     * <p> Two genotypes should have the same number of sites in the same order
     * @param gt1
     * @param gt2
     * @return Return null if the merging is not successful
     */
    public static GenoGrid mergeGenotypesByTaxon(GenoGrid gt1, GenoGrid gt2) {
        if (gt1.getSiteNumber() != gt2.getSiteNumber()) return null;
        int taxaCount = gt1.getTaxaNumber()+gt2.getTaxaNumber();
        BitSet[][] bArray = new BitSet[taxaCount][3];
        String[] taxa = new String[taxaCount];
        int cnt = 0;
        for (int i = 0; i < gt1.getTaxaNumber(); i++) {
            for (int j = 0; j < gt1.genoTaxon[0].length; j++) {
                bArray[cnt][j] = gt1.genoTaxon[i][j];
            }
            taxa[cnt] = gt1.getTaxonName(i);
            cnt++;
        }
        for (int i = 0; i < gt2.getTaxaNumber(); i++) {
            for (int j = 0; j < gt2.genoTaxon[0].length; j++) {
                bArray[cnt][j] = gt2.genoTaxon[i][j];
            }
            taxa[cnt] = gt2.getTaxonName(i);
            cnt++;
        }
        BiSNP[] nsnps = new BiSNP[gt1.getSiteNumber()];
        for (int i = 0; i < nsnps.length; i++) {
            nsnps[i] = gt1.snps[i].replicateWithoutFeature();
        }
        GenoGrid gt = new GenoGrid(bArray, GenoGrid.Direction.ByTaxon, taxa, nsnps);
        return gt;
    }

    /**
     * Convert {@link GenoGrid} to {@link GenoRows}
     * @param gt
     * @return
     */
    public static GenoRows getConvertedGenotype(GenoGrid gt) {
        SiteGeno[] geno = new SiteGeno[gt.getSiteNumber()];
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            short chr = gt.getChromosome(i);
            int pos = gt.getPosition(i);
            char refBase = gt.getReferenceAlleleBase(i);
            char altBase = gt.getAlternativeAlleleBase(i);
            BitSet phase1 = gt.genoSite[i][0];
            BitSet phase2 = gt.genoSite[i][1];
            BitSet missingP = gt.genoSite[i][2];
            geno[i] = new SiteGeno(chr, pos, refBase, altBase, null, phase1, phase2, missingP, gt.getTaxaNumber());
        }
        return new GenoRows(geno, gt.taxa);
    }

    /**
     * Convert {@link GenoGrid} to {@link GenoRows}
     * @param gt
     * @return
     */
    public static GenoGrid getConvertedGenotype(GenoRows gt) {
        BitSet[][] bArray = new BitSet[gt.getSiteNumber()][3];
        int cnt = 0;
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            bArray[cnt][0] = gt.geno[i].phase1;
            bArray[cnt][1] = gt.geno[i].phase2;
            bArray[cnt][2] = gt.geno[i].missing;
        }
        BiSNP[] nsnps = new BiSNP[gt.getSiteNumber()];
        for (int i = 0; i < nsnps.length; i++) {
            nsnps[i] = gt.geno[i];
        }
        return new GenoGrid(bArray, GenoGrid.Direction.ByTaxon, gt.taxa, nsnps);
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param gt
     * @param siteIndices
     * @return
     */
    public static GenoRows getSubsetGenotypeBySite(GenoRows gt, int[] siteIndices) {
        SiteGeno[] ge = new SiteGeno[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            ge[i] = gt.geno[siteIndices[i]];
        }
        GenoRows gb = new GenoRows(ge, gt.taxa);
        return gb;
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param gt
     * @param siteIndices
     * @return
     */
    public static GenoGrid getSubsetGenotypeBySite(GenoGrid gt, int[] siteIndices) {
        BitSet[][] bArray = new BitSet[siteIndices.length][3];
        for (int i = 0; i < siteIndices.length; i++) {
            for (int j = 0; j < bArray[0].length; j++) {
                bArray[i][j] = gt.genoSite[siteIndices[i]][j];
            }
        }
        BiSNP[] nsnps = new BiSNP[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            nsnps[i] = gt.snps[siteIndices[i]];
        }
        return new GenoGrid(bArray, GenoGrid.Direction.BySite, gt.taxa, nsnps);
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param gt
     * @param taxaIndices
     * @return
     */
    public static GenoGrid getSubsetGenotypeByTaxon(GenoGrid gt, int[] taxaIndices) {
        BitSet[][] bArray = new BitSet[taxaIndices.length][3];
        for (int i = 0; i < taxaIndices.length; i++) {
            for (int j = 0; j < bArray[0].length; j++) {
                bArray[i][j] = gt.genoTaxon[taxaIndices[i]][j];
            }
        }
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            nTaxa[i] = gt.getTaxonName(taxaIndices[i]);
        }
        BiSNP[] nsnps = new BiSNP[gt.getSiteNumber()];
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            nsnps[i] = gt.snps[i].replicateWithoutFeature();
        }
        return new GenoGrid(bArray, GenoGrid.Direction.ByTaxon, nTaxa, nsnps);
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param gt
     * @param siteIndices
     */
    public static void subsetsGenotypeBySite(GenoRows gt, int[] siteIndices) {
        gt = getSubsetGenotypeBySite(gt, siteIndices);
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param gt
     * @param siteIndices
     */
    public static void subsetsGenotypeBySite(GenoGrid gt, int[] siteIndices) {
        gt = getSubsetGenotypeBySite(gt, siteIndices);
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param gt
     * @param taxaIndices
     */
    public static void subsetsGenotypeByTaxon(GenoGrid gt, int[] taxaIndices) {
        gt = getSubsetGenotypeByTaxon(gt, taxaIndices);
    }
}
