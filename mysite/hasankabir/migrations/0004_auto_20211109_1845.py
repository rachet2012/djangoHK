# Generated by Django 3.2.8 on 2021-11-09 15:45

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('hasankabir', '0003_remove_database_h'),
    ]

    operations = [
        migrations.AddField(
            model_name='database',
            name='bob',
            field=models.FloatField(default=1.5),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='gamma_gas',
            field=models.FloatField(default=0.7),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='gamma_oil',
            field=models.FloatField(default=0.8),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='gamma_wat',
            field=models.FloatField(default=1),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='muob',
            field=models.FloatField(default=0.5),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='pb',
            field=models.FloatField(default=98),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='rsb',
            field=models.FloatField(default=50),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='database',
            name='t_res',
            field=models.FloatField(default=92),
            preserve_default=False,
        ),
    ]