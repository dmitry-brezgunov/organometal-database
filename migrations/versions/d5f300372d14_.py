"""empty message

Revision ID: d5f300372d14
Revises: 130b81d2390f
Create Date: 2025-07-19 16:28:27.417720

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = 'd5f300372d14'
down_revision: Union[str, None] = '130b81d2390f'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('journals', schema=None) as batch_op:
        batch_op.alter_column('short_name',
               existing_type=sa.VARCHAR(),
               nullable=True)
        batch_op.alter_column('abbreviation',
               existing_type=sa.VARCHAR(),
               nullable=True)
        batch_op.alter_column('old_name',
               existing_type=sa.VARCHAR(),
               nullable=True)

    # ### end Alembic commands ###


def downgrade() -> None:
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('journals', schema=None) as batch_op:
        batch_op.alter_column('old_name',
               existing_type=sa.VARCHAR(),
               nullable=False)
        batch_op.alter_column('abbreviation',
               existing_type=sa.VARCHAR(),
               nullable=False)
        batch_op.alter_column('short_name',
               existing_type=sa.VARCHAR(),
               nullable=False)

    # ### end Alembic commands ###
