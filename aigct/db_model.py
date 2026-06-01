from sqlalchemy import Column, String, Float, Integer, ForeignKey
from sqlalchemy.orm import DeclarativeBase, relationship


class Base(DeclarativeBase):
    pass


class VariantEffectSource(Base):
    __tablename__ = "variant_effect_source"

    code = Column(String(20), primary_key=True)
    name = Column(String(100), nullable=False)
    source_type = Column(String(20))
    description = Column(String(500))

    task_aucs = relationship("VariantEffectTaskAuc", back_populates="source")
    gene_aucs = relationship("VariantEffectGeneAuc", back_populates="source")


class VariantEffectTaskAuc(Base):
    __tablename__ = "variant_effect_task_auc"

    task_code = Column(String(15), primary_key=True)
    score_source = Column(
        String(20), ForeignKey("variant_effect_source.code"), primary_key=True
    )
    auc = Column(Float)
    num_positive = Column(Integer)
    num_negative = Column(Integer)

    source = relationship("VariantEffectSource", back_populates="task_aucs")


class VariantEffectGeneAuc(Base):
    __tablename__ = "variant_effect_gene_auc"

    task_code = Column(String(15), primary_key=True)
    score_source = Column(
        String(20), ForeignKey("variant_effect_source.code"), primary_key=True
    )
    gene_symbol = Column(String(20), primary_key=True)
    auc = Column(Float)
    num_positive = Column(Integer)
    num_negative = Column(Integer)

    source = relationship("VariantEffectSource", back_populates="gene_aucs")
