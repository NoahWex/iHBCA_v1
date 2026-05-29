#!/usr/bin/env python3
"""
Step 5: Build Harmonized Donor Metadata (Phase A.2)

Applies mappings from derived_conditions.yaml to create clean donor-level metadata.

Input: unified_donor_metadata.csv (287 donors × 78 columns)
Output: harmonized_donor_metadata.csv (287 donors × harmonized conditions)

Architecture: Orthogonal dimensions (brca_genotype, cancer_history, tissue_indication)
             derived composite (risk_status_binary)

# Harmonization: making data agree that it was always supposed to be friends.
"""

import pandas as pd
import numpy as np
import yaml
import re
from pathlib import Path
from typing import Optional, Dict


# Kumar per-donor dissociation minutes — loaded from dissociation_protocols.yaml.
# Canonical source: iHBCA_publication/publication/assembly/v1/harmonization/dissociation_protocols.yaml
# Not in the iHBCAv1 uploaded object; curated from Kumar et al. 2023 cell-level metadata
# (suspension_dissociation_time, per-donor mode). Three protocols by study design.
_DISSOC_PROTOCOLS_PATH = Path("/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/assembly/v1/harmonization/dissociation_protocols.yaml")

def _load_dissoc_protocols():
    with open(_DISSOC_PROTOCOLS_PATH) as f:
        protocols = yaml.safe_load(f)
    studies = protocols['studies']

    # Per-donor map for Kumar (variable by design)
    kumar_per_donor = {k: v for k, v in studies['kumar']['per_donor'].items()}

    # Study-level constants for consistent-protocol studies
    # Key: ihbca_dataset value; value: minutes (None = no dissociation)
    study_constants = {}
    for study_key, study_data in studies.items():
        if study_key == 'kumar':
            continue  # handled per-donor
        if study_key == 'reed':
            continue  # handled from per-donor column in unified metadata
        minutes = (
            study_data.get('minutes_total') or
            study_data.get('minutes_typical') or
            study_data.get('minutes_collagenase')
        )
        study_constants[study_key] = int(minutes) if minutes is not None else None

    return kumar_per_donor, study_constants

_KUMAR_DISSOC_MINUTES, _STUDY_DISSOC_CONSTANTS = _load_dissoc_protocols()


# =============================================================================
# ORTHOGONAL RISK DIMENSION FUNCTIONS (Phase 1)
# =============================================================================

def extract_brca_genotype(row: pd.Series) -> Optional[str]:
    """
    Extract BRCA genotype from study-specific columns.

    Returns: 'BRCA1', 'BRCA2', 'RAD51C', 'negative', or 'unknown'
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        # Gray has explicit genotype column with mutation details
        val = row.get('gray_Genotype (mutation)')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == 'WT':
                return 'negative'
            elif val_str.startswith('BRCA1'):
                return 'BRCA1'
            elif val_str.startswith('BRCA2'):
                return 'BRCA2'
            elif val_str.startswith('RAD51C'):
                return 'RAD51C'
        return 'unknown'

    elif study == 'nee':
        # Nee has BRCA1 Mutation column with various formats
        val = row.get('nee_BRCA1 Mutation')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.startswith('No (confirmed'):
                return 'negative'
            elif val_str.startswith('BRCA1') or val_str.startswith('BRCA ('):
                # Note: BRCA (c.xxxx) without "1" or "2" - context indicates BRCA1
                # since column is "BRCA1 Mutation"
                return 'BRCA1'
            elif val_str.lower() == 'unknown':
                # Check Patient No. as fallback
                patient_no = row.get('nee_Patient No.')
                if pd.notna(patient_no) and str(patient_no).startswith('BRCA1'):
                    return 'BRCA1'
                return 'unknown'
        return 'unknown'

    elif study.startswith('pal'):
        # Pal has brca1_mutation column (presence = BRCA1)
        # Also check ihbca_risk_status for additional info
        brca_mut = row.get('pal_brca1_mutation')
        risk_status = row.get('ihbca_risk_status')

        if pd.notna(brca_mut) and str(brca_mut).strip() not in ['', 'NaN', 'nan']:
            return 'BRCA1'
        elif pd.notna(risk_status):
            rs = str(risk_status).strip()
            if rs == 'AR':
                return 'unknown'  # Normal tissue, not tested
            elif 'BR1' in rs:
                return 'BRCA1'
            elif 'BR2' in rs:
                return 'BRCA2'
        return 'unknown'

    elif study == 'reed':
        # Reed has explicit brca_status
        val = row.get('reed_brca_status')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == 'BRCA1':
                return 'BRCA1'
            elif val_str == 'BRCA2':
                return 'BRCA2'
            elif val_str in ('WT', 'assume_WT'):
                return 'negative'
            elif val_str == 'unknown':
                return 'unknown'
        return 'unknown'

    elif study == 'kumar':
        # Kumar: no genetic testing data available
        return 'unknown'

    elif study == 'murrow':
        # Murrow: all reduction mammoplasty, no genetic testing
        return 'unknown'

    elif study == 'twigger':
        # Twigger: lactating donors, no genetic testing
        return 'unknown'

    return 'unknown'


def extract_cancer_history(row: pd.Series) -> Optional[str]:
    """
    Extract personal breast cancer history.

    Returns: 'yes', 'no', or 'unknown'
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        # Gray has explicit breast cancer history column
        val = row.get('gray_Breast cancer history')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str in ['no', 'none', '-', '']:
                return 'no'
            elif val_str not in ['', 'nan']:
                return 'yes'
        # Also check ihbca_risk_status for c prefix (cancer)
        risk_status = row.get('ihbca_risk_status')
        if pd.notna(risk_status):
            rs = str(risk_status).strip()
            if rs.startswith('HR-c'):
                return 'yes'
            elif rs == 'AR':
                return 'no'
        return 'unknown'

    elif study == 'nee':
        # Nee has Cancer History column with detailed info
        val = row.get('nee_Cancer History')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str in ['-', '', 'nan', 'none']:
                return 'no'
            elif 'dcis' in val_str or 'invasive' in val_str or 'carcinoma' in val_str:
                return 'yes'
            elif 'family history' in val_str:
                # Family history is not personal cancer history
                return 'no'
            elif 'ovarian' in val_str:
                # Ovarian cancer is not breast cancer
                return 'no'
        return 'unknown'

    elif study.startswith('pal'):
        # Pal: infer from ihbca_risk_status (c prefix = cancer)
        risk_status = row.get('ihbca_risk_status')
        if pd.notna(risk_status):
            rs = str(risk_status).strip()
            if rs.startswith('HR-c'):
                return 'yes'
            elif rs == 'AR':
                return 'no'
            elif rs.startswith('HR-'):
                # HR without c = high risk but no cancer
                return 'no'
        return 'unknown'

    elif study == 'reed':
        # Reed has explicit Previous_breast_cancer_diagnosis
        val = row.get('reed_Previous_breast_cancer_diagnosis')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == 'No':
                return 'no'
            elif val_str == 'Yes':
                return 'yes'
        return 'unknown'

    elif study == 'kumar':
        # Kumar: Cancer Mastectomy = cancer history
        val = row.get('kumar_Tissue_Source')
        if pd.notna(val):
            val_str = str(val).strip()
            if 'Cancer Mastectomy' in val_str:
                return 'yes'
            elif 'Reduction' in val_str or 'Prophylactic' in val_str:
                return 'no'
        return 'unknown'

    elif study == 'murrow':
        # Murrow: all reduction mammoplasty (general population)
        return 'no'

    elif study == 'twigger':
        # Twigger: lactating donors (general population)
        return 'no'

    return 'unknown'


def extract_tissue_indication(row: pd.Series) -> Optional[str]:
    """
    Extract reason for surgery / tissue source.

    Returns: 'reduction', 'prophylactic', 'cancer_mastectomy', 'contralateral', or 'unknown'
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Surgery')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'reductive' in val_str or 'reduction' in val_str:
                return 'reduction'
            elif 'contralateral' in val_str:
                return 'contralateral'
            elif 'prophylactic' in val_str:
                return 'prophylactic'
        return 'unknown'

    elif study == 'nee':
        val = row.get('nee_Tissue Source')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'reduction' in val_str:
                return 'reduction'
            elif 'prophylactic' in val_str or 'prophylatic' in val_str:  # Note typo in data
                return 'prophylactic'
            elif 'contralateral' in val_str:
                return 'contralateral'
        return 'unknown'

    elif study.startswith('pal'):
        # Pal: infer from ihbca_risk_status and surgical notes
        risk_status = row.get('ihbca_risk_status')
        if pd.notna(risk_status):
            rs = str(risk_status).strip()
            if rs == 'AR':
                return 'reduction'
            elif rs.startswith('HR-c'):
                # Cancer history -> likely contralateral or cancer mastectomy
                return 'contralateral'
            elif rs.startswith('HR-'):
                return 'prophylactic'
        return 'unknown'

    elif study == 'reed':
        val = row.get('reed_surgery_type')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'reduction' in val_str:
                return 'reduction'
            elif 'prophylactic' in val_str:
                return 'prophylactic'
            elif 'contralateral' in val_str:
                return 'contralateral'
        return 'unknown'

    elif study == 'kumar':
        val = row.get('kumar_Tissue_Source')
        if pd.notna(val):
            val_str = str(val).strip()
            if 'Reduction' in val_str:
                return 'reduction'
            elif 'Cancer Mastectomy' in val_str:
                return 'contralateral'  # Contralateral tissue from cancer patient
            elif 'Prophylactic' in val_str:
                return 'prophylactic'
        return 'unknown'

    elif study == 'murrow':
        # Murrow: all reduction mammoplasty
        return 'reduction'

    elif study == 'twigger':
        # Twigger donors are lactating women providing milk for cell collection
        # (no surgery indication). Tagged distinctly from surgical-tissue donors.
        return 'milk-derived'

    return 'unknown'


def extract_parity_count(row: pd.Series) -> Optional[int]:
    """
    Extract numeric parity count from study-specific columns.

    Returns: integer parity (0, 1, 2, ...) or None
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Births')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            try:
                return int(float(val))
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'murrow':
        val = row.get('murrow_Parity')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.lower() == 'unknown':
                return None
            try:
                return int(val_str)
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'nee':
        # GXPX format - use existing parser
        val = row.get('nee_Parity')
        return parse_gxpx_parity(val)

    elif study.startswith('pal'):
        # Pal only has binary (Nulliparous/Parous)
        val = row.get('pal_parity')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'nulliparous' in val_str:
                return 0
            # Parous but unknown count - return None for count
            # (binary will still be derived correctly)
        return None

    elif study == 'twigger':
        val = row.get('twigger_Parity')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            try:
                return int(float(val))
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'kumar':
        val = row.get('kumar_Parity')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.lower() == 'unknown':
                return None
            try:
                return int(val_str)
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'reed':
        val = row.get('reed_parity')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.lower() == 'unknown':
                return None
            try:
                return int(val_str)
            except (ValueError, TypeError):
                return None
        return None

    return None


def derive_parity_binary(row: pd.Series) -> Optional[int]:
    """
    Derive binary parity from parity count.

    Returns: 0 (nulliparous), 1 (parous), or None.
    """
    count = row.get('parity_count')

    if pd.notna(count):
        return 0 if count == 0 else 1

    # Fallback for Pal (has binary but not count)
    study = row.get('study')
    if study and str(study).startswith('pal'):
        val = row.get('pal_parity') if 'pal_parity' in row.index else None
        if val is None:
            # Try from original df passed via closure - not available here
            # Use ihbca data instead
            pass

    return None


def extract_age_at_first_birth(row: pd.Series) -> Optional[int]:
    """
    Extract age at first birth where available.

    Returns: integer age or None
    """
    study = row['ihbca_dataset']

    if study == 'nee':
        val = row.get('nee_Age of first birth')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.lower() in ['-', 'unknown', 'nan', '']:
                return None
            try:
                return int(float(val_str))
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'reed':
        # Reed has comma-separated list of ages at each pregnancy
        # Extract first value as age at first birth
        val = row.get('reed_Approximate_age_at_each_last_pregnancy')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.lower() in ['nan', '']:
                return None
            # Split by comma or slash and take first numeric value
            import re
            parts = re.split(r'[,/]', val_str)
            for part in parts:
                part = part.strip()
                # Extract first number from part (handle "27298" typo, "25/26" etc)
                match = re.match(r'^(\d{2})(?:\d{3})?', part)  # Match 2 digits, ignore extra
                if match:
                    age = int(match.group(1))
                    if 10 <= age <= 50:  # Sanity check
                        return age
        return None

    return None


# =============================================================================
# PHASE 4: ETHNICITY
# =============================================================================

def extract_ethnicity_verbatim(row: pd.Series) -> Optional[str]:
    """
    Extract verbatim ethnicity from study-specific columns.

    Returns: study-specific ethnicity string or None
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Race/ethnicity')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            return str(val).strip()
        return None

    elif study == 'nee':
        val = row.get('nee_Ancestry')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            return str(val).strip()
        return None

    elif study == 'kumar':
        val = row.get('kumar_Ethnicity')
        if pd.notna(val) and str(val).strip().lower() not in ['', 'nan', 'unknown']:
            return str(val).strip()
        return None

    elif study == 'reed':
        val = row.get('reed_Ethnicity')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            return str(val).strip()
        return None

    # murrow, pal, twigger: no ethnicity data
    return None


def derive_ethnicity_grouped(row: pd.Series) -> Optional[str]:
    """
    Derive grouped ethnicity from verbatim.

    Returns: 'white', 'black', 'hispanic', 'asian', 'jewish', 'other', or None
    """
    verbatim = row.get('ethnicity_verbatim')
    if pd.isna(verbatim) or verbatim is None:
        return None

    v = str(verbatim).lower()

    # Jewish (check first - may overlap with white)
    if 'jewish' in v or 'ashkenazi' in v:
        return 'jewish'

    # White/Caucasian (check BEFORE asian to avoid "caucasian" -> "asian" substring match)
    if 'caucasian' in v or 'white' in v or 'british' in v or 'irish' in v or 'european' in v:
        return 'white'

    # Hispanic
    if 'hispanic' in v or 'mexican' in v:
        return 'hispanic'

    # Asian (check after caucasian to avoid substring match)
    if 'asian' in v or 'filipino' in v or 'bangladeshi' in v or 'japanese' in v:
        return 'asian'

    # Indian - check separately (could be South Asian or Native American context)
    # In this dataset context, "indian" typically means South Asian
    if 'indian' in v:
        return 'asian'

    # Black
    if 'black' in v or 'african' in v:
        return 'black'

    # Other
    if 'other' in v or 'not_stated' in v:
        return 'other'

    return 'other'  # Default for unrecognized values


# =============================================================================
# PHASE 5: BMI
# =============================================================================

def extract_bmi_continuous(row: pd.Series) -> Optional[float]:
    """
    Extract continuous BMI from study-specific columns.

    Returns: float BMI or None
    """
    study = row['ihbca_dataset']

    if study == 'murrow':
        val = row.get('murrow_BMI')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == 'unknown':
                return None
            try:
                return float(val_str)
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'nee':
        val = row.get('nee_BMI')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str in ['-', '', 'nan', 'NA']:
                return None
            try:
                return float(val_str)
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'kumar':
        val = row.get('kumar_BMI')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == 'unknown':
                return None
            try:
                return float(val_str)
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'reed':
        val = row.get('reed_Body_mass_index')
        if pd.notna(val):
            try:
                return float(val)
            except (ValueError, TypeError):
                return None
        return None

    # gray, pal, twigger: no BMI data
    return None


def derive_bmi_category(row: pd.Series) -> Optional[str]:
    """
    Derive BMI category from continuous value.

    WHO categories:
    - underweight: <18.5
    - normal: 18.5-24.9
    - overweight: 25-29.9
    - obese: >=30
    """
    bmi = row.get('bmi_continuous')
    if pd.isna(bmi) or bmi is None:
        return None

    if bmi < 18.5:
        return 'underweight'
    elif bmi < 25:
        return 'normal'
    elif bmi < 30:
        return 'overweight'
    else:
        return 'obese'


# =============================================================================
# PHASE 6: SAMPLE PRESERVATION
# =============================================================================

def extract_sample_preservation(row: pd.Series) -> Optional[str]:
    """
    Extract sample preservation method.

    Returns: 'fresh', 'frozen', or None
    """
    study = row['ihbca_dataset']

    if study == 'twigger':
        val = row.get('twigger_Fresh_or_Frozen')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'fresh' in val_str:
                return 'fresh'
            elif 'frozen' in val_str:
                return 'frozen'
        return None

    # Most other studies are fresh tissue (default for scRNA-seq)
    # But we don't have explicit documentation, so return None
    return None


def extract_sample_type(row: pd.Series) -> Optional[str]:
    """
    Extract sample type from ihbca_sample_type.

    Returns: 'mixed', 'supernatant', 'organoid', 'organoid_lp', or None
    """
    val = row.get('ihbca_sample_type')
    if pd.notna(val):
        val_str = str(val).strip().lower()
        if val_str == 'mixed':
            return 'mixed'
        elif val_str == 'supernatant':
            return 'supernatant'
        elif val_str == 'organoid lp':
            return 'organoid_lp'
        elif val_str == 'organoid':
            return 'organoid'
    return None


def extract_facs_status(row: pd.Series) -> Optional[str]:
    """
    Extract FACS sorting status from ihbca_FACS_status.

    Returns: 'live_sorted', 'not_sorted', 'cell_type_sorted', 'no_sort', or None
    """
    val = row.get('ihbca_FACS_status')
    if pd.notna(val):
        val_str = str(val).strip().lower()
        if val_str in ('live_sorted', 'not_sorted', 'cell_type_sorted', 'no_sort'):
            return val_str
    return None


def extract_dissociation_minutes(row: pd.Series) -> Optional[int]:
    """
    Extract dissociation/handling time in minutes for all 7 iHBCA studies.

    Reed: per-donor from reed_dissociation_minutes column (5/7/10 min enzymatic).
    Kumar: per-donor from dissociation_protocols.yaml; three protocols by design.
           8 UCI donors unknown (null).
    Nee/Gray/Murrow/Pal: study-level constant from dissociation_protocols.yaml.
    Twigger: null (milk cells; no tissue dissociation).

    Source: publication/assembly/v1/harmonization/dissociation_protocols.yaml

    Returns: integer minutes or None
    """
    study = row['ihbca_dataset']

    if study == 'reed':
        val = row.get('reed_dissociation_minutes')
        if pd.notna(val):
            try:
                return int(float(val))
            except (ValueError, TypeError):
                return None
        return None

    if study == 'kumar':
        donor_id = row.get('ihbca_donor_id')
        if donor_id and donor_id in _KUMAR_DISSOC_MINUTES:
            return _KUMAR_DISSOC_MINUTES[donor_id]
        return None

    # Consistent-protocol studies: return study-level constant
    if study in _STUDY_DISSOC_CONSTANTS:
        return _STUDY_DISSOC_CONSTANTS[study]

    return None


def extract_menopausal_detailed(row: pd.Series) -> Optional[str]:
    """
    Extract detailed menopausal status from study-specific columns.

    Returns: 'pre', 'peri', 'post_natural', 'post_surgical', or None
    """
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Menopause')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'premenopausal' in val_str:
                return 'pre'
            elif 'surgical' in val_str:
                return 'post_surgical'
            elif 'menopause' in val_str:
                return 'post_natural'
        return None

    elif study == 'nee':
        val = row.get('nee_Menstrual Status')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == '-':
                return None
            elif 'premenopause' in val_str:
                return 'pre'
            elif 'menopausal' in val_str:
                # "Menopausal - Age XX" indicates natural menopause
                return 'post_natural'
        return None

    elif study.startswith('pal'):
        val = row.get('pal_menopausal_status')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'pre' in val_str:
                return 'pre'
            elif 'post' in val_str:
                # Pal doesn't distinguish surgical vs natural
                return 'post_natural'
        return None

    elif study == 'kumar':
        val = row.get('kumar_Menopause')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == 'pre':
                return 'pre'
            elif val_str == 'post':
                # Kumar doesn't distinguish surgical vs natural
                return 'post_natural'
            elif val_str == 'unknown':
                return None
        return None

    elif study == 'reed':
        val = row.get('reed_Menopause_status')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == 'Pre':
                return 'pre'
            elif val_str == 'Post_(surgically_induced)':
                return 'post_surgical'
            elif val_str == 'Post':
                return 'post_natural'
            elif val_str == 'Peri/Post':
                # Conservative: assign to peri
                return 'peri'
        return None

    elif study == 'murrow':
        # Murrow: all donors age 19-47 (premenopausal)
        return 'pre'

    elif study == 'twigger':
        # Lactating donors are by definition pre-menopausal.
        return 'pre'

    return None


def derive_menopausal_binary(row: pd.Series) -> Optional[int]:
    """
    Derive binary menopausal status from detailed status.

    Mapping:
    - pre, peri               -> 0
    - post_natural, post_surgical -> 1
    The readable form (pre/peri/post_natural/post_surgical) lives in
    menopausal_status_detailed; this column is the numeric covariate encoding.
    """
    detailed = row.get('menopausal_status_detailed')

    if detailed in ('pre', 'peri'):
        return 0
    elif detailed in ('post_natural', 'post_surgical'):
        return 1
    return None


def derive_risk_genotype_only(row: pd.Series) -> Optional[str]:
    """
    Derive binary risk based on BRCA genotype only.

    Excludes cancer history and tissue indication from risk assignment.
    Use for cleaner genetic predisposition analysis.

    Returns:
    - HR if brca_genotype IN (BRCA1, BRCA2, RAD51C)
    - AR if brca_genotype = negative (confirmed non-carrier)
    - None if brca_genotype = unknown
    """
    brca = row.get('brca_genotype')

    if brca in ('BRCA1', 'BRCA2', 'RAD51C'):
        return 'HR'
    elif brca == 'negative':
        return 'AR'
    else:
        return None


def derive_risk_status_binary(row: pd.Series) -> Optional[str]:
    """
    Derive binary risk status from orthogonal dimensions.

    Rules:
    - HR if: brca_genotype IN (BRCA1, BRCA2, RAD51C) OR
             cancer_history = yes OR
             tissue_indication IN (prophylactic, cancer_mastectomy, contralateral)
    - AR if: tissue_indication = reduction AND no HR flags
    - AR if: general population (no known risk factors) - includes lactating donors
    - None otherwise (truly unknown)
    """
    brca = row.get('brca_genotype')
    cancer = row.get('cancer_history')
    tissue = row.get('tissue_indication')
    study = row.get('study')

    # Any known BRCA mutation = HR
    if brca in ('BRCA1', 'BRCA2', 'RAD51C'):
        return 'HR'

    # Cancer history = HR
    if cancer == 'yes':
        return 'HR'

    # Prophylactic or cancer-related surgery = HR
    if tissue in ('prophylactic', 'cancer_mastectomy', 'contralateral'):
        return 'HR'

    # Reduction mammoplasty with no known risk factors = AR
    if tissue == 'reduction':
        return 'AR'

    # Tested negative with no other risk factors = AR
    if brca == 'negative' and cancer == 'no':
        return 'AR'

    # Twigger: lactating donors from general population = AR
    # (no genetic testing, no surgery - milk collection from healthy donors)
    if study == 'twigger' and cancer == 'no':
        return 'AR'

    # Unknown everything = None
    return None


def build_metadata_notes(row: pd.Series) -> Optional[str]:
    """
    Build catch-all notes for metadata provenance.

    Records why intermediate fields are None when converted from 'unknown'.
    General-purpose annotation column — extend for any metadata ambiguity.
    """
    # Metadata provenance: because "trust me" isn't a valid citation
    notes = []
    study = row['study']

    if row.get('brca_genotype') == 'unknown':
        if study in ('kumar', 'murrow', 'twigger', 'pal'):
            notes.append("brca: not_tested")
        else:
            notes.append("brca: status_unknown")

    if row.get('cancer_history') == 'unknown':
        notes.append("cancer_hx: not_recorded")

    if study == 'twigger':
        # Twigger donors contribute milk-derived cells (no surgery indication).
        # Recorded explicitly so the note survives even when tissue_indication
        # carries the dedicated 'milk-derived' value rather than 'unknown'.
        notes.append("tissue: milk_collection")
    elif row.get('tissue_indication') == 'unknown':
        notes.append("tissue: not_recorded")

    return '; '.join(notes) if notes else None


def parse_gxpx_parity(value: str) -> Optional[int]:
    """
    Parse GXPX notation to extract parity (P value).

    Examples:
        G0P0M0A0 -> 0
        G3P3M0A0 -> 3
        G2P1 -> 1
        - -> None
    """
    if pd.isna(value) or str(value).strip() in ['-', '', 'nan', 'NA']:
        return None
    match = re.search(r'P(\d+)', str(value))
    if match:
        return int(match.group(1))
    return None


def harmonize_age_continuous(row: pd.Series) -> Optional[float]:
    """Derive continuous age from study-specific columns."""
    study = row['ihbca_dataset']

    # Column mapping by study
    col_map = {
        'gray': 'gray_Age',
        'murrow': 'murrow_Age',
        'nee': 'nee_Age',
        'twigger': 'twigger_Maternal_age',
        'pal': 'pal_patient_age',
        'reed': 'ihbca_age',
    }

    # Kumar excluded (binary only)
    if study == 'kumar':
        return None

    col = col_map.get(study)
    if col and col in row.index:
        val = row[col]
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA', '-']:
            try:
                return float(val)
            except (ValueError, TypeError):
                return None
    return None


def harmonize_age_binary(row: pd.Series, cutoff: int = 50) -> Optional[str]:
    """Derive binary age (young/old) from study-specific columns."""
    study = row['ihbca_dataset']

    # Kumar has native binary encoding
    if study == 'kumar':
        val = row.get('kumar_Age')
        if pd.notna(val):
            val_str = str(val).strip().upper()
            if val_str == 'Y':
                return 'young'
            elif val_str == 'O':
                return 'old'
        return None

    # Other studies: derive from continuous
    age = harmonize_age_continuous(row)
    if age is not None:
        return 'young' if age < cutoff else 'old'
    return None


def harmonize_parity_binary(row: pd.Series) -> Optional[int]:
    """Derive binary parity (0=nulliparous, 1=parous) from study-specific columns."""
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Births')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            try:
                births = float(val)
                return 0 if births == 0 else 1
            except (ValueError, TypeError):
                return None

    elif study == 'murrow':
        val = row.get('murrow_Parity')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == '0':
                return 0
            elif val_str in ['1', '2', '3', '4']:
                return 1
            elif val_str == 'unknown':
                return None
        return None

    elif study == 'nee':
        val = row.get('nee_Parity')
        parity = parse_gxpx_parity(val)
        if parity is not None:
            return 0 if parity == 0 else 1
        return None

    elif study == 'twigger':
        val = row.get('twigger_Parity')
        if pd.notna(val) and str(val).strip() not in ['', 'nan', 'NA']:
            try:
                parity = float(val)
                return 0 if parity == 0 else 1
            except (ValueError, TypeError):
                return None
        return None

    elif study == 'pal':
        val = row.get('pal_parity')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'nulliparous' in val_str:
                return 0
            elif 'parous' in val_str:
                return 1
        return None

    elif study == 'kumar':
        val = row.get('kumar_Parity')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == '0':
                return 0
            elif val_str == '1':
                return 1
            elif val_str.lower() == 'unknown':
                return None
        return None

    elif study == 'reed':
        val = row.get('ihbca_parous')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str in ['true', '1', '1.0']:
                return 1
            elif val_str in ['false', '0', '0.0']:
                return 0
        return None

    return None


def harmonize_risk_status(row: pd.Series) -> Optional[str]:
    """Derive binary risk status (AR/HR) from study-specific columns."""
    study = row['ihbca_dataset']

    # CRITICAL: Nee uses native column (iHBCA labels are WRONG)
    if study == 'nee':
        val = row.get('nee_BRCA1 Mutation')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str.startswith('No (confirmed'):
                return 'AR'
            elif val_str.startswith('BRCA'):
                return 'HR'
            elif val_str.lower() == 'unknown':
                # "Unknown" means unknown specific mutation, not unknown BRCA status
                # Fall back to Patient No. which encodes BRCA status
                patient_no = row.get('nee_Patient No.')
                if pd.notna(patient_no) and str(patient_no).startswith('BRCA1'):
                    return 'HR'
                return None
        return None

    # All other studies use ihbca_risk_status
    val = row.get('ihbca_risk_status')
    if pd.notna(val):
        val_str = str(val).strip()
        if val_str == 'AR':
            return 'AR'
        elif val_str.startswith('HR'):
            return 'HR'
    return None


def harmonize_menopausal_status(row: pd.Series) -> Optional[int]:
    """Derive binary menopausal status (0=pre, 1=post) from study-specific columns."""
    study = row['ihbca_dataset']

    if study == 'gray':
        val = row.get('gray_Menopause')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'premenopausal' in val_str:
                return 0
            elif 'menopause' in val_str or 'surgical' in val_str:
                return 1
        return None

    elif study == 'nee':
        val = row.get('nee_Menstrual Status')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == '-':
                return None
            elif 'premenopause' in val_str:
                return 0
            elif 'menopausal' in val_str:
                return 1
        return None

    elif study == 'pal':
        val = row.get('pal_menopausal_status')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if 'pre' in val_str:
                return 0
            elif 'post' in val_str:
                return 1
        return None

    elif study == 'kumar':
        val = row.get('kumar_Menopause')
        if pd.notna(val):
            val_str = str(val).strip().lower()
            if val_str == 'pre':
                return 0
            elif val_str == 'post':
                return 1
            elif val_str == 'unknown':
                return None
        return None

    elif study == 'reed':
        val = row.get('reed_Menopause_status')
        if pd.notna(val):
            val_str = str(val).strip()
            if val_str == 'Pre':
                return 0
            elif val_str in ['Post', 'Peri/Post', 'Post_(surgically_induced)']:
                return 1
        return None

    elif study == 'murrow':
        return 0

    return None


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Build harmonized donor metadata")
    parser.add_argument("--input", default=None,
                        help="Path to unified_donor_metadata.csv. Defaults to ../outputs/unified_metadata/unified_donor_metadata.csv relative to this script.")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory. Defaults to ../outputs/harmonized_metadata/ relative to this script.")
    args = parser.parse_args()

    # Paths — explicit args take priority; fall back to legacy relative layout
    script_dir = Path(__file__).parent
    base_dir = script_dir.parent
    input_path = Path(args.input) if args.input else base_dir / 'outputs' / 'unified_metadata' / 'unified_donor_metadata.csv'
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / 'outputs' / 'harmonized_metadata'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load existing harmonized metadata for backward compatibility check
    existing_path = output_dir / 'harmonized_donor_metadata.csv'
    existing_df = None
    if existing_path.exists():
        existing_df = pd.read_csv(existing_path)
        print(f"Loaded existing harmonized metadata for validation: {len(existing_df)} donors")

    # Load data
    print(f"Loading {input_path}")
    df = pd.read_csv(input_path)
    # Strip leading/trailing whitespace from donor IDs. Some upstream entries
    # carry stray spaces (e.g. ' Ctrl_ Pt1') that would otherwise break joins
    # on ihbca_donor_id.
    df['ihbca_donor_id'] = df['ihbca_donor_id'].str.strip()
    print(f"Loaded: {len(df)} donors × {len(df.columns)} columns")

    # Apply harmonization
    print("Applying harmonization mappings...")

    harmonized = pd.DataFrame({
        'ihbca_donor_id': df['ihbca_donor_id'],
        'study': df['ihbca_dataset'],
    })

    # Age
    print("  - age_continuous")
    harmonized['age_continuous'] = df.apply(harmonize_age_continuous, axis=1)
    print("  - age_binary")
    harmonized['age_binary'] = df.apply(harmonize_age_binary, axis=1)

    # Parity (count + binary + age at first birth)
    print("  - parity_count")
    harmonized['parity_count'] = df.apply(extract_parity_count, axis=1)
    print("  - parity_binary (derived from count, fallback to legacy)")
    # First derive from count
    harmonized['parity_binary'] = harmonized.apply(derive_parity_binary, axis=1)
    # Fallback to legacy for rows where count is None but legacy has value (e.g., Pal)
    legacy_parity = df.apply(harmonize_parity_binary, axis=1)
    mask = harmonized['parity_binary'].isna() & legacy_parity.notna()
    harmonized.loc[mask, 'parity_binary'] = legacy_parity[mask]
    print("  - age_at_first_birth")
    harmonized['age_at_first_birth'] = df.apply(extract_age_at_first_birth, axis=1)

    # === ORTHOGONAL RISK DIMENSIONS (Phase 1) ===
    print("\n  === Orthogonal Risk Dimensions ===")
    print("  - brca_genotype")
    harmonized['brca_genotype'] = df.apply(extract_brca_genotype, axis=1)
    print("  - cancer_history")
    harmonized['cancer_history'] = df.apply(extract_cancer_history, axis=1)
    print("  - tissue_indication")
    harmonized['tissue_indication'] = df.apply(extract_tissue_indication, axis=1)

    # === PROVENANCE NOTES (capture reasons before converting 'unknown' → None) ===
    print("  - metadata_notes (provenance annotations)")
    harmonized['metadata_notes'] = harmonized.apply(build_metadata_notes, axis=1)
    n_annotated = harmonized['metadata_notes'].notna().sum()
    print(f"    {n_annotated}/{len(harmonized)} donors have provenance notes")

    # Convert intermediate 'unknown' strings to None for robustness
    # Derivation functions handle None identically to 'unknown' — verified in provenance audit
    for col in ['brca_genotype', 'cancer_history', 'tissue_indication']:
        n_unknown = (harmonized[col] == 'unknown').sum()
        harmonized[col] = harmonized[col].replace('unknown', None)
        if n_unknown > 0:
            print(f"    {col}: converted {n_unknown} 'unknown' → None")

    # Derived risk status (from orthogonal dimensions)
    print("  - risk_status_binary (derived from dimensions)")
    harmonized['risk_status_binary'] = harmonized.apply(derive_risk_status_binary, axis=1)

    # Genotype-only risk (for cleaner genetic analysis)
    print("  - risk_genotype_only (BRCA carriers vs confirmed non-carriers)")
    harmonized['risk_genotype_only'] = harmonized.apply(derive_risk_genotype_only, axis=1)

    # Legacy risk status (for comparison)
    print("  - risk_status_binary_legacy (original logic)")
    harmonized['risk_status_binary_legacy'] = df.apply(harmonize_risk_status, axis=1)

    # Menopausal status (detailed + binary)
    print("  - menopausal_status_detailed (pre/peri/post_natural/post_surgical)")
    harmonized['menopausal_status_detailed'] = df.apply(extract_menopausal_detailed, axis=1)
    print("  - menopausal_status_binary (derived from detailed)")
    harmonized['menopausal_status_binary'] = harmonized.apply(derive_menopausal_binary, axis=1)

    # === PHASE 4: ETHNICITY ===
    print("\n  === Ethnicity ===")
    print("  - ethnicity_verbatim")
    harmonized['ethnicity_verbatim'] = df.apply(extract_ethnicity_verbatim, axis=1)
    print("  - ethnicity_grouped (derived)")
    harmonized['ethnicity_grouped'] = harmonized.apply(derive_ethnicity_grouped, axis=1)

    # === PHASE 5: BMI ===
    print("\n  === BMI ===")
    print("  - bmi_continuous")
    harmonized['bmi_continuous'] = df.apply(extract_bmi_continuous, axis=1)
    print("  - bmi_category (derived)")
    harmonized['bmi_category'] = harmonized.apply(derive_bmi_category, axis=1)

    # === PHASE 6: SAMPLE METADATA ===
    print("\n  === Sample Metadata ===")
    print("  - sample_preservation")
    harmonized['sample_preservation'] = df.apply(extract_sample_preservation, axis=1)
    print("  - sample_type")
    harmonized['sample_type'] = df.apply(extract_sample_type, axis=1)
    print("  - facs_status")
    harmonized['facs_status'] = df.apply(extract_facs_status, axis=1)
    print("  - dissociation_minutes")
    harmonized['dissociation_minutes'] = df.apply(extract_dissociation_minutes, axis=1)

    # =========================================================================
    # VALIDATION SECTION
    # =========================================================================
    print("\n" + "=" * 70)
    print("VALIDATION")
    print("=" * 70)

    # 1. Backward compatibility check: risk_status_binary
    print("\n--- BACKWARD COMPATIBILITY: risk_status_binary ---")
    mismatch_mask = harmonized['risk_status_binary'] != harmonized['risk_status_binary_legacy']
    mismatches = harmonized[mismatch_mask].copy()

    if len(mismatches) == 0:
        print("✓ risk_status_binary: 100% match with legacy logic")
    else:
        print(f"⚠ risk_status_binary: {len(mismatches)} mismatches with legacy logic")
        print("\nMismatch details:")
        for _, row in mismatches.iterrows():
            print(f"  {row['ihbca_donor_id']} ({row['study']}): "
                  f"new={row['risk_status_binary']} vs legacy={row['risk_status_binary_legacy']}")
            print(f"    brca={row['brca_genotype']}, cancer={row['cancer_history']}, tissue={row['tissue_indication']}")

    # 2. Existing file comparison (if available)
    if existing_df is not None:
        print("\n--- COMPARISON WITH EXISTING FILE ---")
        # Merge on ihbca_donor_id
        compare = harmonized.merge(
            existing_df[['ihbca_donor_id', 'risk_status_binary']],
            on='ihbca_donor_id',
            suffixes=('_new', '_existing'),
            how='left'
        )
        exist_mismatch = compare[compare['risk_status_binary_new'] != compare['risk_status_binary_existing']]
        if len(exist_mismatch) == 0:
            print("✓ risk_status_binary: 100% match with existing file")
        else:
            print(f"⚠ risk_status_binary: {len(exist_mismatch)} mismatches with existing file")
            for _, row in exist_mismatch.iterrows():
                print(f"  {row['ihbca_donor_id']}: new={row['risk_status_binary_new']} vs existing={row['risk_status_binary_existing']}")

    # 3. Kumar-specific validation (primary test case per plan)
    print("\n--- KUMAR VALIDATION ---")
    kumar = harmonized[harmonized['study'] == 'kumar'].copy()
    kumar_unified = df[df['ihbca_dataset'] == 'kumar'][['ihbca_donor_id', 'kumar_Tissue_Source']].copy()
    kumar = kumar.merge(kumar_unified, on='ihbca_donor_id')

    print("Kumar tissue_indication mapping:")
    kumar_cross = pd.crosstab(kumar['kumar_Tissue_Source'], kumar['tissue_indication'])
    print(kumar_cross.to_string())

    print("\nKumar risk derivation:")
    kumar_risk = pd.crosstab(
        [kumar['kumar_Tissue_Source']],
        [kumar['brca_genotype'], kumar['cancer_history'], kumar['tissue_indication'], kumar['risk_status_binary']],
        dropna=False
    )
    print(kumar_risk.to_string())

    # 4. Orthogonal dimensions coverage
    print("\n--- ORTHOGONAL DIMENSION COVERAGE ---")
    for col in ['brca_genotype', 'cancer_history', 'tissue_indication']:
        print(f"\n{col}:")
        vc = harmonized[col].value_counts(dropna=False)
        for val, count in vc.items():
            pct = 100 * count / len(harmonized)
            print(f"  {val}: {count} ({pct:.1f}%)")

        print("  By study:")
        for study in harmonized['study'].unique():
            study_data = harmonized[harmonized['study'] == study][col]
            dist = study_data.value_counts().to_dict()
            dist_str = ', '.join(f"{k}={v}" for k, v in dist.items() if pd.notna(k))
            print(f"    {study}: {dist_str}")

    # 5. Genotype-only risk validation
    print("\n--- GENOTYPE-ONLY RISK ---")
    geno_risk = harmonized['risk_genotype_only'].value_counts(dropna=False)
    print("Distribution:")
    for val, count in geno_risk.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")
    print("\nBy study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['risk_genotype_only']
        dist = study_data.value_counts().to_dict()
        n_valid = study_data.notna().sum()
        dist_str = ', '.join(f"{k}={v}" for k, v in dist.items() if pd.notna(k))
        print(f"  {study}: {n_valid}/{len(study_data)} ({dist_str})")

    # 6. Parity count validation
    print("\n--- PARITY COUNT ---")
    parity_counts = harmonized['parity_count'].value_counts(dropna=False).sort_index()
    print("Distribution:")
    for val, count in parity_counts.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")
    print("\nBy study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['parity_count']
        n_valid = study_data.notna().sum()
        if n_valid > 0:
            dist = study_data.value_counts().sort_index().to_dict()
            dist_str = ', '.join(f"{int(k) if pd.notna(k) else 'None'}={v}" for k, v in dist.items() if pd.notna(k))
            print(f"  {study}: {n_valid}/{len(study_data)} ({dist_str})")
        else:
            print(f"  {study}: {n_valid}/{len(study_data)} (binary only)")

    # 7. Age at first birth validation
    print("\n--- AGE AT FIRST BIRTH ---")
    afb = harmonized['age_at_first_birth']
    n_valid = afb.notna().sum()
    print(f"Coverage: {n_valid}/{len(harmonized)} ({100*n_valid/len(harmonized):.1f}%)")
    if n_valid > 0:
        print(f"Range: {int(afb.min())}-{int(afb.max())}")
    print("By study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['age_at_first_birth']
        n_study_valid = study_data.notna().sum()
        if n_study_valid > 0:
            print(f"  {study}: {n_study_valid}/{len(study_data)} (range: {int(study_data.min())}-{int(study_data.max())})")
        else:
            print(f"  {study}: {n_study_valid}/{len(study_data)}")

    # 8. Menopausal status detailed validation
    print("\n--- MENOPAUSAL STATUS DETAILED ---")
    meno_detailed = harmonized['menopausal_status_detailed'].value_counts(dropna=False)
    print("Distribution:")
    for val, count in meno_detailed.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")
    print("\nBy study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['menopausal_status_detailed']
        dist = study_data.value_counts().to_dict()
        n_valid = study_data.notna().sum()
        dist_str = ', '.join(f"{k}={v}" for k, v in dist.items() if pd.notna(k))
        print(f"  {study}: {n_valid}/{len(study_data)} ({dist_str})")

    # 9. Ethnicity validation
    print("\n--- ETHNICITY ---")
    print("Verbatim coverage:")
    n_verb = harmonized['ethnicity_verbatim'].notna().sum()
    print(f"  {n_verb}/{len(harmonized)} ({100*n_verb/len(harmonized):.1f}%)")
    print("\nGrouped distribution:")
    eth_grouped = harmonized['ethnicity_grouped'].value_counts(dropna=False)
    for val, count in eth_grouped.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")
    print("\nBy study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['ethnicity_grouped']
        n_valid = study_data.notna().sum()
        if n_valid > 0:
            dist = study_data.value_counts().to_dict()
            dist_str = ', '.join(f"{k}={v}" for k, v in dist.items() if pd.notna(k))
            print(f"  {study}: {n_valid}/{len(study_data)} ({dist_str})")
        else:
            print(f"  {study}: {n_valid}/{len(study_data)}")

    # 10. BMI validation
    print("\n--- BMI ---")
    bmi = harmonized['bmi_continuous']
    n_bmi = bmi.notna().sum()
    print(f"Continuous coverage: {n_bmi}/{len(harmonized)} ({100*n_bmi/len(harmonized):.1f}%)")
    if n_bmi > 0:
        print(f"Range: {bmi.min():.1f} - {bmi.max():.1f}")
    print("\nCategory distribution:")
    bmi_cat = harmonized['bmi_category'].value_counts(dropna=False)
    for val, count in bmi_cat.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")
    print("\nBy study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['bmi_continuous']
        n_valid = study_data.notna().sum()
        if n_valid > 0:
            print(f"  {study}: {n_valid}/{len(study_data)} (range: {study_data.min():.1f}-{study_data.max():.1f})")
        else:
            print(f"  {study}: {n_valid}/{len(study_data)}")

    # 11. Sample metadata validation
    print("\n--- SAMPLE METADATA ---")

    print("\nsample_preservation:")
    samp = harmonized['sample_preservation'].value_counts(dropna=False)
    for val, count in samp.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")

    print("\nsample_type:")
    stype = harmonized['sample_type'].value_counts(dropna=False)
    for val, count in stype.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")

    print("\nfacs_status:")
    facs = harmonized['facs_status'].value_counts(dropna=False)
    for val, count in facs.items():
        pct = 100 * count / len(harmonized)
        print(f"  {val}: {count} ({pct:.1f}%)")

    print("\ndissociation_minutes (Reed + Kumar):")
    dissoc = harmonized['dissociation_minutes']
    n_dissoc = dissoc.notna().sum()
    print(f"  Coverage: {n_dissoc}/{len(harmonized)} ({100*n_dissoc/len(harmonized):.1f}%)")
    if n_dissoc > 0:
        dissoc_dist = dissoc.value_counts().sort_index().to_dict()
        dist_str = ', '.join(f"{int(k)}min={v}" for k, v in dissoc_dist.items() if pd.notna(k))
        print(f"  Values: {dist_str}")
    print("  By study:")
    for study in harmonized['study'].unique():
        study_data = harmonized[harmonized['study'] == study]['dissociation_minutes']
        n_study_valid = study_data.notna().sum()
        print(f"    {study}: {n_study_valid}/{len(study_data)}")

    # 12. Standard coverage validation (all columns)
    print("\n--- STANDARD COVERAGE ---")
    for col in ['age_continuous', 'age_binary', 'parity_count', 'parity_binary', 'age_at_first_birth',
                'risk_status_binary', 'risk_genotype_only', 'menopausal_status_detailed', 'menopausal_status_binary',
                'ethnicity_grouped', 'bmi_continuous', 'bmi_category',
                'sample_preservation', 'sample_type', 'facs_status',
                'dissociation_minutes']:
        n_valid = harmonized[col].notna().sum()
        n_total = len(harmonized)
        print(f"\n{col}:")
        print(f"  Coverage: {n_valid}/{n_total} ({100*n_valid/n_total:.1f}%)")

        if col != 'age_continuous':
            value_counts = harmonized[col].value_counts(dropna=False)
            for val, count in value_counts.items():
                print(f"    {val}: {count}")

        print("  By study:")
        for study in harmonized['study'].unique():
            study_data = harmonized[harmonized['study'] == study][col]
            n_study_valid = study_data.notna().sum()
            n_study_total = len(study_data)
            if col == 'age_continuous':
                if n_study_valid > 0:
                    print(f"    {study}: {n_study_valid}/{n_study_total} (range: {study_data.min():.0f}-{study_data.max():.0f})")
                else:
                    print(f"    {study}: {n_study_valid}/{n_study_total}")
            else:
                dist = study_data.value_counts().to_dict()
                dist_str = ', '.join(f"{k}={v}" for k, v in dist.items() if pd.notna(k))
                print(f"    {study}: {n_study_valid}/{n_study_total} ({dist_str})")

    # =========================================================================
    # OUTPUT
    # =========================================================================

    # Drop legacy column before saving
    harmonized_output = harmonized.drop(columns=['risk_status_binary_legacy'])

    # Reorder columns for clarity
    column_order = [
        'ihbca_donor_id', 'study',
        # Age
        'age_continuous', 'age_binary',
        # Parity
        'parity_count', 'parity_binary', 'age_at_first_birth',
        # Risk (orthogonal dimensions)
        'brca_genotype', 'cancer_history', 'tissue_indication',
        # Risk (derived)
        'risk_status_binary', 'risk_genotype_only',
        # Provenance
        'metadata_notes',
        # Menopausal
        'menopausal_status_detailed', 'menopausal_status_binary',
        # Ethnicity
        'ethnicity_verbatim', 'ethnicity_grouped',
        # BMI
        'bmi_continuous', 'bmi_category',
        # Sample metadata
        'sample_preservation', 'sample_type', 'facs_status',
        'dissociation_minutes'
    ]
    harmonized_output = harmonized_output[column_order]

    # Cast nullable-integer columns so the CSV serializes 0/1/<empty> rather
    # than 0.0/1.0/<empty> from pandas's default object-with-floats coercion.
    # Limited to columns this build redefines as numeric binaries; existing
    # numeric columns (parity_count, age_at_first_birth) are left in their
    # current float serialization to keep the diff minimal.
    for col in ('parity_binary', 'menopausal_status_binary', 'dissociation_minutes'):
        harmonized_output[col] = pd.array(harmonized_output[col], dtype='Int64')

    # Write output
    output_path = output_dir / 'harmonized_donor_metadata.csv'
    print(f"\nWriting {output_path}")
    harmonized_output.to_csv(output_path, index=False)

    print(f"Done. Harmonized metadata saved with {len(harmonized_output)} donors.")

    # Summary
    print("\n=== SUMMARY ===")
    print(f"Total donors: {len(harmonized_output)}")
    print(f"Columns: {list(harmonized_output.columns)}")


if __name__ == '__main__':
    main()
