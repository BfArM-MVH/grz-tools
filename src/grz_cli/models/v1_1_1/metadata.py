from __future__ import annotations

import json
import logging
import typing
from collections.abc import Generator
from datetime import date
from enum import StrEnum
from importlib.resources import files
from itertools import groupby
from pathlib import Path
from typing import Annotated, Any, Self

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StringConstraints,
    model_validator,
)
from pydantic.alias_generators import to_camel

from grz_cli.file_operations import calculate_sha256  # type: ignore

SCHEMA_URL = "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.1.1/GRZ/grz-schema.json"

log = logging.getLogger(__name__)


class StrictBaseModel(BaseModel):
    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        use_enum_values=True,
        alias_generator=to_camel,
    )


class SubmissionType(StrEnum):
    """
    The options are: 'initial' for first submission, 'followup' is for followup submissions, 'addition' for additional submission, 'correction' for correction
    """

    initial = "initial"
    followup = "followup"
    addition = "addition"
    correction = "correction"


class GenomicStudyType(StrEnum):
    """
    whether additional persons are tested as well
    """

    single = "single"
    duo = "duo"
    trio = "trio"


class GenomicStudySubtype(StrEnum):
    """
    whether tumor and/or germ-line are tested
    """

    tumor_only = "tumor-only"
    tumor_germline = "tumor+germline"
    germline_only = "germline-only"


class CoverageType(StrEnum):
    """
    Health insurance providers
    """

    GKV = "GKV"
    PKV = "PKV"
    BG = "BG"
    SEL = "SEL"
    SOZ = "SOZ"
    GPV = "GPV"
    PPV = "PPV"
    BEI = "BEI"
    SKT = "SKT"
    UNK = "UNK"


class DiseaseType(StrEnum):
    """
    Type of the disease
    """

    oncological = "oncological"
    rare = "rare"
    hereditary = "hereditary"


class Submission(StrictBaseModel):
    submission_date: date
    """
    Date of submission in ISO 8601 format YYYY-MM-DD
    """

    submission_type: SubmissionType
    """
    The options are: 'initial' for first submission, 'followup' is for followup submissions, 'addition' for additional submission, 'correction' for correction
    """

    tan_g: Annotated[str, StringConstraints(pattern=r"^[a-fA-F0-9]{64}$")]
    """
    The VNg of the genomic data of the index patient that will be reimbursed --> a unique 32-length byte code represented in a hex string of length 64.
    """

    local_case_id: str
    """
    A local case identifier for synchronizing locally
    """

    coverage_type: CoverageType
    """
    Health insurance provider
    """

    submitter_id: Annotated[str, StringConstraints(pattern=r"^[0-9]{9}$")]
    """
    Institutional ID of the submitter according to §293 SGB V.
    """

    genomic_data_center_id: Annotated[
        str,
        StringConstraints(pattern=r"^(GRZ)[A-Z0-9]{3}[0-9]{3}$"),
    ]
    """
    ID of the genomic data center in the format GRZXXXnnn.
    """

    clinical_data_node_id: Annotated[str, StringConstraints(pattern=r"^(KDK)[A-Z0-9]{3}[0-9]{3}$")]
    """
    ID of the clinical data node in the format KDKXXXnnn.
    """

    disease_type: DiseaseType
    """
    Type of the disease
    """

    genomic_study_type: GenomicStudyType
    """
    whether additional persons are tested as well
    """

    genomic_study_subtype: GenomicStudySubtype
    """
    whether tumor and/or germ-line are tested
    """

    lab_name: str
    """
    Name of the sequencing lab.
    """


class Gender(StrEnum):
    """
    Gender of the donor.
    """

    male = "male"
    female = "female"
    other = "other"
    unknown = "unknown"


class Relation(StrEnum):
    """
    Relationship of the donor in respect to the index patient, e.g. 'index', 'brother', 'mother', etc.
    """

    mother = "mother"
    father = "father"
    brother = "brother"
    sister = "sister"
    child = "child"
    index_ = "index"
    other = "other"


class MvConsentScopeType(StrEnum):
    """
    Consent or refusal to participate and consent, must be indicated for each option listed in the scope of consent.
    """

    permit = "permit"
    deny = "deny"


class MvConsentScopeDomain(StrEnum):
    """
    Scope of consent or revocation.
    """

    mv_sequencing = "mvSequencing"
    re_identification = "reIdentification"
    case_identification = "caseIdentification"


class MvConsentScope(StrictBaseModel):
    """
    The scope of the Modellvorhaben consent given by the donor.
    """

    type_: Annotated[MvConsentScopeType, Field(alias="type")]
    """
    Consent or refusal to participate and consent, must be indicated for each option listed in the scope of consent.
    """

    date: date
    """
    Date of signature of the pilot projects consent; in ISO 8601 format YYYY-MM-DD.
    """

    domain: MvConsentScopeDomain
    """
    Scope of consent or revocation.
    """


class MvConsent(StrictBaseModel):
    presentation_date: date | None = None
    """
    Date of delivery. Date (in ISO 8601 format YYYY-MM-DD) on which the Model Project Declaration of Participation 
    was presented to the patient, unless identical to the date of signature
    """

    version: str
    """
    Version of the declaration of participation. 
    Name and version of the declaration of participation in the MV GenomSeq, e.g.: 
    'Patient Info TE Consent MVGenomSeq vers01'
    """

    scope: list[MvConsentScope]
    """
    Modules of the consent to MV: must have at least a permit of mvSequencing
    """

    @model_validator(mode="after")
    def ensure_mv_sequencing_scope_is_present(self):
        if not any(
            scope.domain == MvConsentScopeDomain.mv_sequencing and scope.type_ == MvConsentScopeType.permit
            for scope in self.scope
        ):
            raise ValueError("Must have at least a permit of mvSequencing")
        return self


class ResearchConsentSchemaVersion(StrEnum):
    v_2025_0_1 = "2025.0.1"


class ResearchConsent(StrictBaseModel):
    """
    Research consents. Multiple declarations of consent are possible! Must be assigned to the respective data sets.
    """

    schema_version: ResearchConsentSchemaVersion
    """
    Schema version of de.medizininformatikinitiative.kerndatensatz.consent
    """

    presentation_date: date | None = None
    """
    Date of the delivery of the research consent in ISO 8601 format (YYYY-MM-DD)
    """

    scope: object
    """
    Scope of the research consent in JSON format following the MII IG Consent v2025 FHIR schema. 
    See 'https://www.medizininformatik-initiative.de/Kerndatensatz/KDS_Consent_V2025/MII-IG-Modul-Consent.html' and 
    'https://packages2.fhir.org/packages/de.medizininformatikinitiative.kerndatensatz.consent'.
    """


class TissueOntology(StrictBaseModel):
    name: str
    """
    Name of the tissue ontology
    """

    version: str
    """
    Version of the tissue ontology
    """


class SampleConservation(StrEnum):
    """
    Sample conservation
    """

    fresh_tissue = "fresh-tissue"
    cryo_frozen = "cryo-frozen"
    ffpe = "ffpe"
    other = "other"
    unknown = "unknown"


class SequenceType(StrEnum):
    """
    Type of sequence (DNA or RNA)
    """

    dna = "dna"
    rna = "rna"


class SequenceSubtype(StrEnum):
    """
    Subtype of sequence (germline, somatic, etc.)
    """

    germline = "germline"
    somatic = "somatic"
    other = "other"
    unknown = "unknown"


class FragmentationMethod(StrEnum):
    """
    Fragmentation method
    """

    sonication = "sonication"
    enzymatic = "enzymatic"
    none = "none"
    other = "other"
    unknown = "unknown"


class LibraryType(StrEnum):
    """
    Library type
    """

    panel = "panel"
    panel_lr = "panel_lr"
    wes = "wes"
    wes_lr = "wes_lr"
    wgs = "wgs"
    wgs_lr = "wgs_lr"
    wxs = "wxs"
    wxs_lr = "wxs_lr"
    other = "other"
    unknown = "unknown"


class EnrichmentKitManufacturer(StrEnum):
    """
    Manufacturer of the enrichment kit
    """

    illumina = "Illumina"
    agilent = "Agilent"
    twist = "Twist"
    neb = "NEB"
    other = "other"
    unknown = "unknown"
    none = "none"


class SequencingLayout(StrEnum):
    """
    The sequencing layout, aka the end type of sequencing.
    """

    single_end = "single-end"
    paired_end = "paired-end"
    reverse = "reverse"
    other = "other"


class TumorCellCountMethod(StrEnum):
    """
    Method used to determine cell count.
    """

    pathology = "pathology"
    bioinformatics = "bioinformatics"
    other = "other"
    unknown = "unknown"


class TumorCellCount(StrictBaseModel):
    """
    Tuple of tumor cell counts and how they were determined.
    """

    count: Annotated[float, Field(ge=0.0, le=100.0)]
    """
    Tumor cell count in %
    """

    method: TumorCellCountMethod
    """
    Method used to determine cell count.
    """


class CallerUsedItem(StrictBaseModel):
    name: str
    """
    Name of the caller used
    """

    version: str
    """
    Version of the caller used
    """


class FileType(StrEnum):
    """
    Type of the file; if BED file is submitted, only 1 file is allowed.
    """

    bam = "bam"
    vcf = "vcf"
    bed = "bed"
    fastq = "fastq"


class ChecksumType(StrEnum):
    """
    Type of checksum algorithm used
    """

    sha256 = "sha256"


class ReadOrder(StrEnum):
    """
    Indicates the read order for paired-end reads.
    """

    r1 = "R1"
    r2 = "R2"


class File(StrictBaseModel):
    file_path: str
    """
    Path relative to the submission root, e.g.: sequencing_data/patient_001/patient_001_dna.bam
    """

    file_type: FileType
    """
    Type of the file; if BED file is submitted, only 1 file is allowed.
    """

    read_length: Annotated[int, Field(strict=True, ge=0)] | None = None
    """
    The read length; in the case of long-read sequencing it is the rounded average read length.
    """

    checksum_type: ChecksumType | None = ChecksumType.sha256
    """
    Type of checksum algorithm used
    """

    file_checksum: Annotated[str, StringConstraints(pattern=r"^[a-fA-F0-9]{64}$")]
    """
    checksum of the file
    """

    file_size_in_bytes: Annotated[int, Field(strict=True, ge=0)]
    """
    Size of the file in bytes
    """

    read_order: ReadOrder | None = None
    """
    Indicates the read order for paired-end reads.
    """

    flowcell_id: str | None = None
    """
    Indicates the flow cell.
    """

    lane_id: str | None = None
    """
    Indicates the lane
    """

    @model_validator(mode="after")
    def ensure_read_length_is_present_for_bam_and_fastq(self):
        if self.file_type in {FileType.bam, FileType.fastq} and self.read_length is None:
            raise ValueError(f"Read length missing for file '{self.file_path}' of type '{self.file_type}'.")
        return self

    def validate_data(self, local_file_path: Path) -> Generator[str]:
        """
        Validates whether the provided file matches this metadata.

        :param local_file_path: Path to the actual file (resolved if symlinked)
        :return: Generator of errors
        """
        # Resolve file path
        local_file_path = local_file_path.resolve()

        # Check if path exists
        if not local_file_path.exists():
            yield f"{str(self.file_path)} does not exist!"
            # Return here as following tests cannot work
            return

        # Check if path is a file
        if not local_file_path.is_file():
            yield f"{str(self.file_path)} is not a file!"
            # Return here as following tests cannot work
            return

        # Check if the checksum is correct
        if self.checksum_type == "sha256":
            calculated_checksum = calculate_sha256(local_file_path)
            if self.file_checksum != calculated_checksum:
                yield (
                    f"{str(self.file_path)}: Checksum mismatch! "
                    f"Expected: '{self.file_checksum}', calculated: '{calculated_checksum}'."
                )
        else:
            yield (
                f"{str(self.file_path)}: Unsupported checksum type: {self.checksum_type}. "
                f"Supported types: {[e.value for e in ChecksumType]}"
            )

        # Check file size
        if self.file_size_in_bytes != local_file_path.stat().st_size:
            yield (
                f"{str(self.file_path)}: File size mismatch! "
                f"Expected: '{self.file_size_in_bytes}', observed: '{local_file_path.stat().st_size}'."
            )

    def encrypted_file_path(self):
        return self.file_path + ".c4gh"


class PercentBasesAboveQualityThreshold(StrictBaseModel):
    """Percentage of bases with a specified minimum quality threshold"""

    minimum_quality: Annotated[float, Field(strict=True, ge=0.0)]
    """The minimum quality score threshold"""

    percent: Annotated[float, Field(strict=True, ge=0.0, le=100.0)]
    """
    Percentage of bases with a specified minimum quality threshold, according to https://www.bfarm.de/SharedDocs/Downloads/DE/Forschung/modellvorhaben-genomsequenzierung/Qs-durch-GRZ.pdf?__blob=publicationFile
    """


class ReferenceGenome(StrEnum):
    """
    Reference genome used according to the Genome Reference Consortium (https://www.ncbi.nlm.nih.gov/grc)
    """

    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"


class SequenceData(StrictBaseModel):
    bioinformatics_pipeline_name: str
    """
    Name of the bioinformatics pipeline used
    """

    bioinformatics_pipeline_version: str
    """
    Version or commit hash of the bioinformatics pipeline
    """

    reference_genome: ReferenceGenome
    """
    Reference genome used according to the Genome Reference Consortium (https://www.ncbi.nlm.nih.gov/grc)
    """

    percent_bases_above_quality_threshold: PercentBasesAboveQualityThreshold
    """
    Percentage of bases with a specified minimum quality threshold
    """

    mean_depth_of_coverage: Annotated[float, Field(strict=True, ge=0.0)]
    """
    Mean depth of coverage
    """

    min_coverage: Annotated[float, Field(strict=True, ge=0.0)]
    """
    Minimum coverage
    """

    targeted_regions_above_min_coverage: Annotated[float, Field(strict=True, ge=0.0, le=1.0)]
    """
    Fraction of targeted regions that are above minimum coverage
    """

    non_coding_variants: bool
    """
    The analysis includes non-coding variants -> true or false
    """

    caller_used: list[CallerUsedItem]
    """
    Caller that is used in the pipeline
    """

    files: list[File]
    """
    List of files generated and required in this analysis.
    """

    def contains_files(self, file_type: FileType) -> bool:
        return any(f.file_type == file_type for f in self.files)

    def list_files(self, file_type: FileType) -> list[File]:
        return [f for f in self.files if f.file_type == file_type]


class LabDatum(StrictBaseModel):
    lab_data_name: str
    """
    Name/ID of the biospecimen e.g. 'Blut DNA normal'
    """

    tissue_ontology: TissueOntology

    tissue_type_id: str
    """
    Tissue ID according to the ontology in use.
    """

    tissue_type_name: str
    """
    Tissue name according to the ontology in use.
    """

    sample_date: date
    """
    Date of sample in ISO 8601 format YYYY-MM-DD
    """

    sample_conservation: SampleConservation
    """
    Sample conservation
    """

    sequence_type: SequenceType
    """
    Type of sequence (DNA or RNA)
    """

    sequence_subtype: SequenceSubtype
    """
    Subtype of sequence (germline, somatic, etc.)
    """

    fragmentation_method: FragmentationMethod
    """
    Fragmentation method
    """

    library_type: LibraryType
    """
    Library type
    """

    library_prep_kit: str
    """
    Name/version of the library prepkit
    """

    library_prep_kit_manufacturer: str
    """
    Library prep kit manufacturer
    """

    sequencer_model: str
    """
    Name/version of the sequencer model
    """

    sequencer_manufacturer: str
    """
    Sequencer manufacturer
    """

    kit_name: str
    """
    Name/version of the sequencing kit
    """

    kit_manufacturer: str
    """
    Sequencing kit manufacturer
    """

    enrichment_kit_manufacturer: EnrichmentKitManufacturer
    """
    Manufacturer of the enrichment kit
    """

    enrichment_kit_description: str
    """
    Name/version of the enrichment kit
    """

    barcode: str
    """
    The barcode used or 'na'
    """

    sequencing_layout: SequencingLayout
    """
    The sequencing layout, aka the end type of sequencing.
    """

    tumor_cell_count: list[TumorCellCount] | None = None
    """
    Tuple of tumor cell counts and how they were determined.
    """

    sequence_data: SequenceData | None = None
    """
    Sequence data generated from the wet lab experiment.
    """

    def has_sequence_data(self) -> bool:
        return self.sequence_data is not None

    @model_validator(mode="after")
    def validate_sequencing_setup(self) -> Self:
        if self.library_type in {LibraryType.wxs, LibraryType.wxs_lr} and self.sequence_type != SequenceType.rna:
            raise ValueError(
                f"Error in lab datum '{self.lab_data_name}': "
                f"WXS requires RNA sequencing, but got '{self.sequence_type}'."
            )
        return self


class Donor(StrictBaseModel):
    donor_pseudonym: str
    """
    A unique identifier given by the Leistungserbringer for each donor of a single, duo or trio sequencing; 
    the donorPseudonym needs to be identifiable by the Leistungserbringer 
    in case of changes to the consents by one of the donors. 
    For Index patient use TanG.
    """

    gender: Gender
    """
    Gender of the donor.
    """

    relation: Relation
    """
    Relationship of the donor in respect to the index patient, e.g. 'index', 'brother', 'mother', etc.
    """

    mv_consent: MvConsent

    research_consents: list[ResearchConsent]
    """
    Research consents. Multiple declarations of consent are possible! Must be assigned to the respective data sets.
    """

    lab_data: list[LabDatum]
    """
    Lab data related to the donor.
    """

    @model_validator(mode="after")
    def warn_empty_sequence_data(self):
        for lab_datum in self.lab_data:
            if not lab_datum.has_sequence_data():
                log.warning(
                    f"No sequence data found for lab datum '{lab_datum.lab_data_name}' in donor '{self.donor_pseudonym}'. "
                    "Is this a submission without sequence data?"
                )
        return self

    @model_validator(mode="after")
    def validate_target_bed_files_exist(self):
        """
        Check if the submission has the required bed files for panel sequencing.
        """
        lib_types = {
            LibraryType.panel,
            LibraryType.wes,
            LibraryType.wxs,
            LibraryType.panel_lr,
            LibraryType.wes_lr,
            LibraryType.wxs_lr,
        }

        for lab_datum in self.lab_data:
            if (
                lab_datum.has_sequence_data()
                and lab_datum.library_type in lib_types
                and not lab_datum.sequence_data.contains_files(FileType.bed)
            ):
                raise ValueError(
                    f"BED file missing for lab datum '{lab_datum.lab_data_name}' in donor '{self.donor_pseudonym}'."
                )

        return self

    @model_validator(mode="after")
    def validate_vcf_file_exists(self):
        """
        Check if there is a VCF file
        """
        for lab_datum in self.lab_data:
            if lab_datum.has_sequence_data() and not lab_datum.sequence_data.contains_files(FileType.vcf):
                raise ValueError(
                    f"VCF file missing for lab datum '{lab_datum.lab_data_name}' in donor '{self.donor_pseudonym}'."
                )

        return self

    @model_validator(mode="after")
    def validate_fastq_file_exists(self):  # noqa: C901
        """
        Check if there is a FASTQ file
        """
        for lab_datum in self.lab_data:
            if not lab_datum.has_sequence_data():
                # Skip if no sequence data is present
                continue
            fastq_files = lab_datum.sequence_data.list_files(FileType.fastq)

            if len(fastq_files) == 0:
                raise ValueError("No FASTQ file found!")
            elif lab_datum.sequencing_layout == SequencingLayout.paired_end:
                # check if read order is specified
                for i in fastq_files:
                    if i.read_order is None:
                        raise ValueError(
                            f"Error in lab datum '{lab_datum.lab_data_name}' of donor '{self.donor_pseudonym}': "
                            f"No read order specified for FASTQ file '{i.file_path}'!"
                        )

                key = lambda f: (f.flowcell_id, f.lane_id)
                fastq_files.sort(key=key)
                for _key, group in groupby(fastq_files, key):
                    flowcell_id = _key[0]
                    lane_id = _key[1]
                    files = list(group)

                    # separate R1 and R2 files
                    fastq_r1_files = [f for f in files if f.read_order == ReadOrder.r1]
                    fastq_r2_files = [f for f in files if f.read_order == ReadOrder.r2]

                    # check that there are exactly one R1 and on R2 file present
                    if len(fastq_r1_files) > 1:
                        raise ValueError(
                            f"Error in lab datum '{lab_datum.lab_data_name}' of donor '{self.donor_pseudonym}': "
                            f"Paired end sequencing layout but multiple R1 files for flowcell id '{flowcell_id}', lane id '{lane_id}'!"
                        )
                    elif len(fastq_r1_files) < 1:
                        raise ValueError(
                            f"Error in lab datum '{lab_datum.lab_data_name}' of donor '{self.donor_pseudonym}': "
                            f"Paired end sequencing layout but missing R1 file for flowcell id '{flowcell_id}', lane id '{lane_id}'!"
                        )

                    if len(fastq_r2_files) > 1:
                        raise ValueError(
                            f"Error in lab datum '{lab_datum.lab_data_name}' of donor '{self.donor_pseudonym}': "
                            f"Paired end sequencing layout but multiple R2 files for flowcell id '{flowcell_id}', lane id '{lane_id}'!"
                        )
                    elif len(fastq_r2_files) < 1:
                        raise ValueError(
                            f"Error in lab datum '{lab_datum.lab_data_name}' of donor '{self.donor_pseudonym}': "
                            f"Paired end sequencing layout but missing R2 file for flowcell id '{flowcell_id}', lane id '{lane_id}'!"
                        )

        return self


class GrzSubmissionMetadata(StrictBaseModel):
    """
    General metadata schema for submissions to the GRZ
    """

    schema_: Annotated[str, Field(alias="$schema")] = SCHEMA_URL

    submission: Submission

    donors: list[Donor]
    """
    List of donors including the index patient.
    """

    @model_validator(mode="after")
    def check_schema(self):
        if self.schema_ != SCHEMA_URL:
            log.warning(f"Unknown GRZ metadata schema URL: {self.schema_}")
        return self

    @model_validator(mode="after")
    def validate_donor_count(self):
        """
        Check whether the submission has the required number of donors based on the study type.
        """
        study_type = self.submission.genomic_study_type

        match study_type:
            case GenomicStudyType.single:
                # Check if the submission has at least one donor
                if not self.donors:
                    raise ValueError("At least one donor is required for a single study.")
            case GenomicStudyType.duo:
                # Check if the submission has at least two donors
                if len(self.donors) < 2:
                    raise ValueError("At least two donors are required for a duo study.")
            case GenomicStudyType.trio:
                # Check if the submission has at least three donors
                if len(self.donors) < 3:
                    raise ValueError("At least three donors are required for a trio study.")

        return self

    @model_validator(mode="after")
    def check_for_tumor_cell_count(self):
        """
        Check if oncology samples have tumor cell counts.
        """
        for donor in self.donors:
            pseudonym = donor.donor_pseudonym
            for lab_datum in donor.lab_data:
                if lab_datum.sequence_subtype == SequenceSubtype.somatic and lab_datum.tumor_cell_count is None:
                    raise ValueError(
                        f"Missing tumor cell count for donor '{pseudonym}', lab datum '{lab_datum.lab_data_name}'!"
                    )

        return self

    @model_validator(mode="after")
    def check_duplicate_lab_data_names(self):
        """
        Check if the submission contains lab data with the same name within one donor.
        """
        for donor in self.donors:
            pseudonym = donor.donor_pseudonym
            lab_data_names = set()
            for lab_datum in donor.lab_data:
                if lab_datum.lab_data_name in lab_data_names:
                    raise ValueError(f"Duplicate lab datum '{lab_datum.lab_data_name}' in donor '{pseudonym}'")
                else:
                    lab_data_names.add(lab_datum.lab_data_name)

        return self

    @model_validator(mode="after")
    def validate_thresholds(self):
        """
        Check if the submission meets the minimum mean coverage requirements.
        """
        threshold_definitions = _load_thresholds()

        for donor in self.donors:
            for lab_datum in donor.lab_data:
                key = (
                    self.submission.genomic_study_subtype,
                    lab_datum.library_type,
                    lab_datum.sequence_subtype,
                )
                thresholds = threshold_definitions.get(key)
                if thresholds is None:
                    allowed_combinations = sorted(list(threshold_definitions.keys()))
                    allowed_combinations = "\n".join([f"  - {combination}" for combination in allowed_combinations])
                    names = (
                        "submission.genomicStudySubtype",
                        "labData.libraryType",
                        "labData.sequenceSubtype",
                    )
                    info = dict(zip(names, key, strict=True))
                    log.warning(
                        f"No thresholds for the specified combination {info} found (donor {donor.donor_pseudonym})!\n"
                        f"Valid combinations:\n{allowed_combinations}.\n"
                        f"See https://www.bfarm.de/SharedDocs/Downloads/DE/Forschung/modellvorhaben-genomsequenzierung/Qs-durch-GRZ.pdf?__blob=publicationFile for more details.\n"
                        f"Skipping threshold validation."
                    )
                    continue

                _check_thresholds(donor, lab_datum, thresholds)

        return self

    @model_validator(mode="after")
    def validate_reference_genome_compatibility(self):
        reference_genomes = {
            (donor.donor_pseudonym, lab_datum.lab_data_name): lab_datum.sequence_data.reference_genome
            for donor in self.donors
            for lab_datum in donor.lab_data
        }
        unique_reference_genomes = set(reference_genomes.values())
        if len(unique_reference_genomes) > 1:
            raise ValueError(
                f"Incompatible reference genomes found: {unique_reference_genomes}.\n"
                f"Reference genomes must be consistent within a submission.\n"
                f"Reference genomes: {reference_genomes}"
            )

        return self


def _check_thresholds(donor: Donor, lab_datum: LabDatum, thresholds: dict[str, Any]):
    if not lab_datum.has_sequence_data():
        # Skip if no sequence data is present; warning issues in the validator `warn_empty_sequence_data` of `Donor`.
        return
    pseudonym = donor.donor_pseudonym
    lab_data_name = lab_datum.lab_data_name
    # mypy cannot reason about the `has_sequence_data` check
    sequence_data = typing.cast(SequenceData, lab_datum.sequence_data)

    mean_depth_of_coverage_t = thresholds.get("meanDepthOfCoverage")
    mean_depth_of_coverage_v = sequence_data.mean_depth_of_coverage
    if mean_depth_of_coverage_t and mean_depth_of_coverage_v < mean_depth_of_coverage_t:
        raise ValueError(
            f"Mean depth of coverage for donor '{pseudonym}', lab datum '{lab_data_name}' "
            f"below threshold: {mean_depth_of_coverage_v} < {mean_depth_of_coverage_t}"
        )

    read_length_t = thresholds.get("readLength")
    for f in sequence_data.list_files(FileType.fastq) + sequence_data.list_files(FileType.bam):
        read_length_v = f.read_length
        if read_length_t and read_length_v < read_length_t:
            raise ValueError(
                f"Read length for donor '{pseudonym}', lab datum '{lab_data_name}' "
                f"below threshold: {read_length_v} < {read_length_t}"
            )

    if t := thresholds.get("targetedRegionsAboveMinCoverage"):
        min_coverage_t = t.get("minCoverage")
        min_coverage_v = sequence_data.min_coverage

        fraction_above_t = t.get("fractionAbove")
        fraction_above_v = sequence_data.targeted_regions_above_min_coverage

        if min_coverage_t and min_coverage_v < min_coverage_t:
            raise ValueError(
                f"Minimum coverage for donor '{pseudonym}', lab datum '{lab_data_name}' "
                f"below threshold: {min_coverage_v} < {min_coverage_t}"
            )
        if fraction_above_t and fraction_above_v < fraction_above_t:
            raise ValueError(
                f"Fraction of targeted regions above minimum coverage for donor '{pseudonym}', "
                f"lab datum '{lab_data_name}' below threshold: "
                f"{fraction_above_v} < {fraction_above_t}"
            )


type Thresholds = dict[tuple[str, str, str], dict[str, Any]]


def _load_thresholds() -> Thresholds:
    threshold_definitions = json.load(
        files("grz_cli").joinpath("resources", "thresholds.json").open("r", encoding="utf-8")
    )
    threshold_definitions = {
        (d["genomicStudySubtype"], d["libraryType"], d["sequenceSubtype"]): d["thresholds"]
        for d in threshold_definitions
    }
    return threshold_definitions
